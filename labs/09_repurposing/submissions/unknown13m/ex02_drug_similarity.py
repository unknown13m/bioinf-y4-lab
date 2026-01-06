from pathlib import Path
import pandas as pd
import itertools

HANDLE = "unknown13m"
BASE = Path(f"labs/09_repurposing/submissions/{HANDLE}")
DRUG_GENE_PATH = Path("data/sample/drug_gene_tp53_subset.csv")


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def main():
    BASE.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DRUG_GENE_PATH)

    # detectez coloanele
    drug_col = next((c for c in df.columns if "drug" in c.lower()), df.columns[0])
    gene_col = next((c for c in df.columns if "gene" in c.lower()), df.columns[1])

    df = df[[drug_col, gene_col]].dropna()
    df.columns = ["drug", "gene"]

    # facem mapping cu drug-set(genes)
    drug2genes = {}
    for drug, gene in df.itertuples(index=False):
        drug2genes.setdefault(str(drug).strip(), set()).add(str(gene).strip())

    drugs = sorted(drug2genes.keys())

    rows = []
    for d1, d2 in itertools.combinations(drugs, 2):
        sim = jaccard(drug2genes[d1], drug2genes[d2])
        rows.append({"drug1": d1, "drug2": d2, "jaccard": sim})

    out = pd.DataFrame(rows).sort_values("jaccard", ascending=False)
    out_path = BASE / f"drug_similarity_{HANDLE}.csv"
    out.to_csv(out_path, index=False)

    print("[OK] Wrote:", out_path)
    print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
