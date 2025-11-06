import json, requests, pandas as pd
from pathlib import Path

CBIO = "https://www.cbioportal.org/api"

def fetch_mrna(gene_symbols, study="brca_tcga_pan_can_atlas_2018", profile="brca_tcga_pan_can_atlas_2018_mrna"):
    samples = requests.get(f"{CBIO}/studies/{study}/samples", headers={"accept":"application/json"}).json()
    sample_ids = [s["sampleId"] for s in samples]
    payload = {
        "entrezGeneIds": [],
        "geneSymbols": gene_symbols,
        "sampleIds": sample_ids
    }
    r = requests.post(f"{CBIO}/molecular-profiles/{profile}/mRNA/fetch", json=payload, headers={"content-type":"application/json","accept":"application/json"})
    r.raise_for_status()
    rows = r.json()
    rec = {}
    for row in rows:
        g = row["geneSymbol"]
        for sid, val in row["values"].items():
            rec.setdefault(g, {})[sid] = val
    return pd.DataFrame.from_dict(rec, orient='index')

if __name__ == "__main__":
    Path("data").mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv("data/mapping_table.csv")
    genes = sorted(set(g for line in mapping["genes"] for g in line.split(";")))
    df = fetch_mrna(genes)
    df.to_csv("data/tcga_brca_expr.csv")
    print("Wrote data/tcga_brca_expr.csv")
