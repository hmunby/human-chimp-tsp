import json, time, urllib.request

BASE = "https://gtexportal.org/api/v2/association/dyneqtl"
GENE = "ENSG00000171195.11"          # MUC7
VARIANTS = [
    ("chr4_70344431_A_G_b38", "rs34735123"),
    ("chr4_70344481_C_T_b38", "rs13101613"),
    ("chr4_70344604_A_T_b38", "rs67708234"),
    ("chr4_70344619_T_C_b38", "rs13134277"),
    ("chr4_70345253_A_G_b38", "rs13148087"),
]
TISSUES = ["Minor_Salivary_Gland", "Skin_Sun_Exposed_Lower_leg"]

def get(url):
    for a in range(4):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            body = e.read().decode()[:300]
            if e.code in (400, 404):
                return {"__http__": e.code, "__body__": body}
            if a == 3: return {"__http__": e.code, "__body__": body}
            time.sleep(2*(a+1))
        except Exception as e:
            if a == 3: return {"__error__": str(e)}
            time.sleep(2*(a+1))

out = []
for tis in TISSUES:
    for vid, rs in VARIANTS:
        u = (f"{BASE}?variantId={vid}&gencodeId={GENE}"
             f"&tissueSiteDetailId={tis}&datasetId=gtex_v10")
        d = get(u)
        rec = {"tissue": tis, "rsid": rs, "variantId": vid}
        if isinstance(d, dict) and ("__http__" in d or "__error__" in d):
            rec["tested"] = False
            rec["reason"] = d.get("__body__") or d.get("__error__")
        else:
            rec["tested"] = True
            for k in ("pValue", "nes", "maf", "error", "tissueSiteDetailId",
                      "pValueThreshold", "geneSymbol"):
                if isinstance(d, dict) and k in d:
                    rec[k] = d[k]
            rec["_raw_keys"] = sorted(d.keys()) if isinstance(d, dict) else None
        out.append(rec)
        time.sleep(0.3)

json.dump(out, open("gtex_dyneqtl_v10.json", "w"), indent=1)
for r in out:
    print(r["tissue"], r["rsid"], "tested=", r["tested"],
          "p=", r.get("pValue"), "nes=", r.get("nes"), "maf=", r.get("maf"),
          r.get("reason", "")[:120] if r.get("reason") else "")
