#!/usr/bin/env python3
"""Build the human MEMBRANE-glycoprotein gene set with coordinates.

There is no dedicated UniProt keyword for "membrane glycoprotein"; it is the intersection
of KW-0325 (Glycoprotein) and KW-0472 (Membrane). We fetch that intersection from the UniProt
REST API (human, reviewed) and then subset the existing all-glycoprotein coords table
(human_glycoproteins_KW0325_with_coords.tsv, built by fetch_glycoprotein_coords.py) by UniProt
accession -- so the coordinates are identical to the parent set, no BioMart re-run needed.

Definitions and human-reviewed counts (for reference):
  Glycoprotein (KW-0325)                     4726
  Glycoprotein AND Membrane (KW-0472)        3537   <- this set
  Glycoprotein AND Cell membrane (KW-1003)   2419
  Glycoprotein AND Transmembrane (KW-0812)   3150

Output: human_membrane_glycoproteins_KW0325_KW0472_with_coords.tsv
        (same columns as the parent coords file; 3374 genes with coords, 166 accessions
         had no ENSG/coord mapping in the parent set)
"""
import requests, pandas as pd
from io import StringIO

GD = "resources/glycoprotein_distribution"
PARENT = f"{GD}/human_glycoproteins_KW0325_with_coords.tsv"
OUT = f"{GD}/human_membrane_glycoproteins_KW0325_KW0472_with_coords.tsv"

print("Fetching KW-0325 & KW-0472 (human, reviewed) accessions from UniProt...")
r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                 params={"query": "keyword:KW-0325 AND keyword:KW-0472 AND organism_id:9606 AND reviewed:true",
                         "fields": "accession,gene_primary", "format": "tsv"}, timeout=120)
r.raise_for_status()
mg = pd.read_csv(StringIO(r.text), sep="\t")
acc = set(mg["Entry"])
print(f"  {len(acc)} membrane-glycoprotein accessions")

coords = pd.read_csv(PARENT, sep="\t")
sub = coords[coords["uniprot_accession"].isin(acc)].copy()
sub.to_csv(OUT, sep="\t", index=False)
missing = acc - set(coords["uniprot_accession"])
print(f"  {sub['hgnc_symbol'].nunique()} genes with coords ({len(sub)} rows) -> {OUT}")
print(f"  {len(missing)} accessions had no coord mapping in the parent set")
