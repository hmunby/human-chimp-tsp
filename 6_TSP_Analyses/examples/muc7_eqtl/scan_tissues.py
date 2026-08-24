import json, time, urllib.request

BASE = "https://www.ebi.ac.uk/eqtl/api/v2"
GENES = {"MUC7": "ENSG00000171195", "CABS1": "ENSG00000145309",
         "SMR3A": "ENSG00000109208", "SMR3B": "ENSG00000171201",
         "HTN1": "ENSG00000126550", "HTN3": "ENSG00000205649"}

def get(url):
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if a == 3: return {"__error__": str(e)}
            time.sleep(2*(a+1))

ds = get(f"{BASE}/datasets?study_label=GTEx&quant_method=ge&size=1000")
out = {}
for d in ds:
    sg, did, n = d["sample_group"], d["dataset_id"], d["sample_size"]
    out[sg] = {"dataset": did, "n": n}
    for nm, gid in GENES.items():
        r = get(f"{BASE}/datasets/{did}/associations?gene_id={gid}&size=1")
        out[sg][nm] = (None if isinstance(r, dict)
                       else r[0].get("median_tpm"))
        time.sleep(0.15)

json.dump(out, open("tissue_scan.json", "w"), indent=1)
print("done", len(out))
