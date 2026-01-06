from pathlib import Path
import pandas as pd
import networkx as nx

HANDLE = "unknown13m"
BASE = Path(f"labs/09_repurposing/submissions/{HANDLE}")

DRUG_GENE_PATH = Path("data/sample/drug_gene_tp53_subset.csv")
DISEASE_GENES_PATH = BASE / f"disease_genes_{HANDLE}.txt"


def main():
    #facem load drug-gene edges
    df = pd.read_csv(DRUG_GENE_PATH)
    drug_col = next((c for c in df.columns if "drug" in c.lower()), df.columns[0])
    gene_col = next((c for c in df.columns if "gene" in c.lower()), df.columns[1])
    df = df[[drug_col, gene_col]].dropna()
    df.columns = ["drug", "gene"]

    # incarc lista de disease genes
    disease_genes = []
    for line in DISEASE_GENES_PATH.read_text(encoding="utf-8").splitlines():
        g = line.strip()
        if g:
            disease_genes.append(g)

    # fac un graf bipartit(drug-gene)
    G = nx.Graph()
    for drug, gene in df.itertuples(index=False):
        drug = str(drug).strip()
        gene = str(gene).strip()
        G.add_node(drug, bipartite="drug")
        G.add_node(gene, bipartite="gene")
        G.add_edge(drug, gene)

    drugs = sorted({d for d in df["drug"].unique()})

    # 4) calculam mean distance de la fiecare drug la disease genes
    #daca un disease gene e deconectat, il punem penalizare cu o distanta mare
    PENALTY = 999

    rows = []
    for drug in drugs:
        dists = []
        for g in disease_genes:
            if drug in G and g in G:
                try:
                    dist = nx.shortest_path_length(G, source=drug, target=g)
                except nx.NetworkXNoPath:
                    dist = PENALTY
            else:
                dist = PENALTY
            dists.append(dist)

        mean_dist = sum(dists) / len(dists) if dists else PENALTY
        rows.append(
            {
                "drug": drug,
                "mean_distance": mean_dist,
                "targets": ",".join(sorted({x for x in G.neighbors(drug)})) if drug in G else "",
            }
        )

    out = pd.DataFrame(rows).sort_values("mean_distance", ascending=True)
    out_path = BASE / f"drug_priority_{HANDLE}.csv"
    out.to_csv(out_path, index=False)

    print("[OK] Wrote:", out_path)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
