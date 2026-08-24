#!/usr/bin/env python3
"""GWAS traits reported at each TSP, for the Table 1 / supplementary-table GWAS column.

`opentargets_assign.py` caches the full OpenTargets variant record for every TSP, which
includes each variant's GWAS credible sets and their study ids -- but not the study TRAITS.
This resolves those study ids to trait names and emits one row per (TSP, study).

    tsp_all_gwas_hits.long.tsv    region, variant_id_hg38, rsid, study_id, gwas_trait

Traits are keyed to the variant's hg38 coordinate, so downstream tables match on position
rather than region id and the file stays valid for any TSP set.

Study ids are cached in ot_study_cache.json as the script runs, so a re-run in the
same directory is offline. The cache is not bundled.

Usage:
    python build_gwas_hits.py --tsps ../2_ClusteringLD/tsps.pooled.hg19.txt
"""
import argparse
import json
import os
import time
import urllib.error
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
VARIANT_CACHE = os.path.join(HERE, "ot_variant_cache.json")
STUDY_CACHE = os.path.join(HERE, "ot_study_cache.json")
BATCH = 25


def post(query, retries=4):
    body = json.dumps({"query": query}).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(OT_API, data=body,
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
            print(f"    retry {attempt + 1} after {type(e).__name__}")
    return {}


def fetch_traits(study_ids, cache):
    """Resolve study ids to trait names, in aliased batches, caching every response."""
    todo = [s for s in study_ids if s not in cache]
    print(f"  studies: {len(study_ids)} total, {len(todo)} not cached")
    for i in range(0, len(todo), BATCH):
        chunk = todo[i:i + BATCH]
        subs = "\n".join(
            f'  s{j}: study(studyId: "{sid}") {{ id traitFromSource '
            f'diseases {{ name }} }}' for j, sid in enumerate(chunk))
        data = post("query {\n" + subs + "\n}") or {}
        d = data.get("data") or {}
        for j, sid in enumerate(chunk):
            rec = d.get(f"s{j}") or {}
            trait = rec.get("traitFromSource")
            if not trait:
                dis = rec.get("diseases") or []
                trait = dis[0]["name"] if dis and dis[0].get("name") else None
            cache[sid] = trait
        print(f"    resolved {min(i + BATCH, len(todo))}/{len(todo)}", flush=True)
        with open(STUDY_CACHE, "w") as fh:
            json.dump(cache, fh)
    return cache


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tsps", default=os.path.join(HERE, "../2_ClusteringLD/tsps.pooled.hg19.txt"))
    ap.add_argument("--variant-cache", default=VARIANT_CACHE)
    ap.add_argument("--out", default=os.path.join(HERE, "tsp_all_gwas_hits.long.tsv"))
    a = ap.parse_args()

    vc = json.load(open(a.variant_cache))
    t = pd.read_csv(a.tsps, sep="\t", low_memory=False)
    if "POSITION_HG38" not in t.columns and "POS_HG38" in t.columns:
        t = t.rename(columns={"POS_HG38": "POSITION_HG38"})
    t = t.dropna(subset=["POSITION_HG38"])

    # collect every GWAS credible set for every TSP
    hits, studies, not_indexed, uncached = [], set(), 0, 0
    for r in t.itertuples():
        vid = f"{r.CHROM_HG38}_{int(r.POSITION_HG38)}_{r.REF}_{r.ALT}"
        if vid not in vc:
            uncached += 1                  # never queried -- would need an API call
            continue
        rec = vc[vid]
        if rec is None:
            not_indexed += 1               # queried, but absent from the OpenTargets index
            continue
        rsid = (rec.get("rsIds") or [None])[0]
        for cs in ((rec.get("credibleSets") or {}).get("rows") or []):
            if cs.get("studyType") != "gwas":
                continue
            sid = (cs.get("study") or {}).get("id")
            if not sid:
                continue
            studies.add(sid)
            hits.append(dict(region=r.region, variant_id_hg38=vid, rsid=rsid, study_id=sid))
    print(f"  TSPs: {len(t)}   in the OpenTargets index: {len(t) - not_indexed - uncached}"
          f"   not indexed: {not_indexed}   never queried: {uncached}")
    if uncached:
        print("    -> run opentargets_assign.py first; those variants have no GWAS data yet")
    print(f"  GWAS credible-set hits: {len(hits)} across {len(studies)} studies")

    cache = json.load(open(STUDY_CACHE)) if os.path.exists(STUDY_CACHE) else {}
    cache = fetch_traits(sorted(studies), cache)

    df = pd.DataFrame(hits)
    df["gwas_trait"] = df["study_id"].map(cache)
    n_unres = int(df["gwas_trait"].isna().sum())
    df = df.dropna(subset=["gwas_trait"]).drop_duplicates()
    df.to_csv(a.out, sep="\t", index=False)
    print(f"  unresolved traits dropped: {n_unres}")
    print(f"  wrote {a.out}: {len(df)} rows, "
          f"{df['variant_id_hg38'].nunique()} TSPs with >=1 trait, "
          f"{df['gwas_trait'].nunique()} distinct traits")


if __name__ == "__main__":
    main()
