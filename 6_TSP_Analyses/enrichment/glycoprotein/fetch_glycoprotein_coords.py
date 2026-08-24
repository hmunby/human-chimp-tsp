"""
Fetch all human glycoproteins (KW-0325) from UniProt and retrieve:
  - ENSG gene IDs (from ENST via Ensembl BioMart)
  - Genomic coordinates in hg38 (GRCh38)
  - Genomic coordinates in hg19 (GRCh37)

Fix: BioMart queries use POST instead of GET to avoid 414 URI Too Large errors.

Output: human_glycoproteins_KW0325_with_coords.tsv
"""

import requests
import pandas as pd
from io import StringIO
import time

# ── 1. Fetch from UniProt ────────────────────────────────────────────────────

print("Step 1: Fetching KW-0325 human glycoproteins from UniProt...")

uniprot_url = "https://rest.uniprot.org/uniprotkb/stream"
uniprot_params = {
    "query": "keyword:KW-0325 AND organism_id:9606 AND reviewed:true",
    "fields": "accession,gene_names,gene_primary,xref_ensembl",
    "format": "tsv"
}

r = requests.get(uniprot_url, params=uniprot_params)
r.raise_for_status()

uniprot_df = pd.read_csv(StringIO(r.text), sep="\t")
print(f"  -> {len(uniprot_df)} UniProt entries retrieved")

# ── 2. Extract ENST IDs ──────────────────────────────────────────────────────

print("\nStep 2: Extracting Ensembl transcript IDs (ENST)...")

enst_col = [c for c in uniprot_df.columns if "ensembl" in c.lower()]
if not enst_col:
    raise ValueError(f"Could not find Ensembl column. Available: {uniprot_df.columns.tolist()}")
enst_col = enst_col[0]

rows = []
for _, row in uniprot_df.iterrows():
    raw = str(row[enst_col])
    if raw.strip() in ("", "nan"):
        continue
    enst_ids = [x.strip().split(".")[0]
                for x in raw.split(";") if x.strip().startswith("ENST")]
    for enst in enst_ids:
        rows.append({
            "uniprot_accession": row.get("Entry", ""),
            "gene_name_primary": row.get("Gene Names (primary)", ""),
            "gene_names_all":    row.get("Gene Names", ""),
            "enst_id":           enst
        })

enst_df = pd.DataFrame(rows).drop_duplicates("enst_id")
all_enst = enst_df["enst_id"].tolist()
print(f"  -> {len(all_enst)} unique ENST IDs extracted")

# ── Helper: BioMart POST query ───────────────────────────────────────────────

def biomart_query_post(host, dataset, attributes, filter_name, values,
                       chunk_size=200, max_retries=5, retry_wait=10):
    """
    Query Ensembl BioMart using POST requests (avoids 414 URI Too Large).
    Queries in chunks and concatenates results.
    Retries each chunk up to max_retries times on transient server errors.
    """
    url = f"{host}/biomart/martservice"
    all_results = []
    n_chunks = (len(values) - 1) // chunk_size + 1

    for i in range(0, len(values), chunk_size):
        chunk = values[i : i + chunk_size]
        value_str = ",".join(chunk)
        chunk_num = i // chunk_size + 1

        attr_xml = "\n    ".join(f'<Attribute name="{a}" />' for a in attributes)
        xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1"
       uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="{dataset}" interface="default">
    <Filter name="{filter_name}" value="{value_str}"/>
    {attr_xml}
  </Dataset>
</Query>"""

        # Retry loop for transient BioMart/DB errors
        for attempt in range(1, max_retries + 1):
            try:
                resp = requests.post(url, data={"query": xml}, timeout=120)
                resp.raise_for_status()

                if resp.text.strip().startswith("Query ERROR"):
                    raise ValueError(f"BioMart query error: {resp.text[:200]}")

                chunk_df = pd.read_csv(StringIO(resp.text), sep="\t")
                all_results.append(chunk_df)
                print(f"  -> chunk {chunk_num}/{n_chunks}: {len(chunk_df)} rows")
                time.sleep(0.3)
                break  # success — exit retry loop

            except Exception as e:
                if attempt < max_retries:
                    print(f"  chunk {chunk_num}/{n_chunks} attempt {attempt} failed "
                          f"({str(e)[:120]}). Retrying in {retry_wait}s...")
                    time.sleep(retry_wait)
                else:
                    print(f"  chunk {chunk_num}/{n_chunks} failed after {max_retries} "
                          f"attempts — skipping. Error: {str(e)[:120]}")

    if not all_results:
        return pd.DataFrame()
    return pd.concat(all_results, ignore_index=True)

STD_CHROMS = [str(c) for c in list(range(1, 23)) + ["X", "Y", "MT"]]

# ── 3. BioMart hg38 (GRCh38) ────────────────────────────────────────────────

print("\nStep 3: Querying Ensembl BioMart (hg38 / GRCh38)...")

hg38_raw = biomart_query_post(
    host="https://www.ensembl.org",
    dataset="hsapiens_gene_ensembl",
    attributes=[
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "hgnc_symbol",
    ],
    filter_name="ensembl_transcript_id",
    values=all_enst,
    chunk_size=200,
)

if not hg38_raw.empty:
    hg38_raw = hg38_raw[hg38_raw["Chromosome/scaffold name"].isin(STD_CHROMS)]
    hg38_df = hg38_raw.rename(columns={
        "Transcript stable ID":     "enst_id",
        "Gene stable ID":           "ensg_id",
        "Chromosome/scaffold name": "chrom_hg38",
        "Gene start (bp)":          "start_hg38",
        "Gene end (bp)":            "end_hg38",
        "Strand":                   "strand_hg38",
        "HGNC symbol":              "hgnc_symbol",
    }).drop_duplicates(subset=["enst_id", "ensg_id"])
    print(f"  -> {hg38_df['ensg_id'].nunique()} unique ENSG IDs in hg38")
else:
    hg38_df = pd.DataFrame()
    print("  Warning: no hg38 results returned")

# ── 4. BioMart hg19 (GRCh37) ────────────────────────────────────────────────

print("\nStep 4: Querying Ensembl BioMart (hg19 / GRCh37)...")

hg19_raw = biomart_query_post(
    host="https://grch37.ensembl.org",
    dataset="hsapiens_gene_ensembl",
    attributes=[
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
    ],
    filter_name="ensembl_transcript_id",
    values=all_enst,
    chunk_size=200,
)

if not hg19_raw.empty:
    hg19_raw = hg19_raw[hg19_raw["Chromosome/scaffold name"].isin(STD_CHROMS)]
    hg19_df = hg19_raw.rename(columns={
        "Transcript stable ID":     "enst_id",
        "Gene stable ID":           "ensg_id_hg19",
        "Chromosome/scaffold name": "chrom_hg19",
        "Gene start (bp)":          "start_hg19",
        "Gene end (bp)":            "end_hg19",
        "Strand":                   "strand_hg19",
    }).drop_duplicates(subset=["enst_id", "ensg_id_hg19"])
    print(f"  -> {hg19_df['ensg_id_hg19'].nunique()} unique ENSG IDs in hg19")
else:
    hg19_df = pd.DataFrame()
    print("  Warning: no hg19 results returned")

# ── 5. Merge & collapse to gene level ───────────────────────────────────────

print("\nStep 5: Merging all data...")

merged = hg38_df.copy() if not hg38_df.empty else enst_df.copy()

# Add UniProt info
merged = merged.merge(
    enst_df[["enst_id", "uniprot_accession", "gene_name_primary", "gene_names_all"]],
    on="enst_id", how="left"
)

# Add hg19 coordinates
if not hg19_df.empty:
    hg19_slim = (hg19_df[["enst_id", "chrom_hg19", "start_hg19", "end_hg19", "strand_hg19"]]
                 .drop_duplicates("enst_id"))
    merged = merged.merge(hg19_slim, on="enst_id", how="left")

# Collapse to one row per ENSG gene
hg19_agg = ({"chrom_hg19": "first", "start_hg19": "first",
             "end_hg19": "first", "strand_hg19": "first"}
            if "chrom_hg19" in merged.columns else {})

gene_level = (
    merged.groupby("ensg_id", as_index=False).agg({
        "hgnc_symbol":       "first",
        "gene_name_primary": "first",
        "gene_names_all":    "first",
        "uniprot_accession": lambda x: ";".join(x.dropna().unique()),
        "enst_id":           lambda x: ";".join(x.dropna().unique()),
        "chrom_hg38":        "first",
        "start_hg38":        "first",
        "end_hg38":          "first",
        "strand_hg38":       "first",
        **hg19_agg
    })
)

# ── 6. Save ──────────────────────────────────────────────────────────────────

out_path = "human_glycoproteins_KW0325_with_coords.tsv"
gene_level.to_csv(out_path, sep="\t", index=False)

print(f"\nDone! Output saved to: {out_path}")
print(f"Total genes: {len(gene_level)}")
print("\nColumn summary:")
for col in gene_level.columns:
    non_null = gene_level[col].notna().sum()
    print(f"  {col:30s} {non_null}/{len(gene_level)} non-null")
