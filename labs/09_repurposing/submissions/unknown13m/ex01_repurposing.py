from pathlib import Path
import pandas as pd
import networkx as nx

HANDLE = "unknown13m"
BASE = Path(f"labs/09_repurposing/submissions/{HANDLE}")

DRUG_GENE_PATH = Path("data/sample/drug_gene_tp53_subset.csv")


def main():
    BASE.mkdir(parents=True, exist_ok=True)

    # incarc drug-gene edges
    df = pd.read_csv(DRUG_GENE_PATH)

    # incerc sa inferez cele 2 coloane (drug, gene)
    cols = [c.lower() for c in df.columns]
    # pattern uri comune: drug/gene, drug_name/gene etc etc
    drug_col = None
    gene_col = None
    for c in df.columns:
        cl = c.lower()
        if drug_col is None and "drug" in cl:
            drug_col = c
        if gene_col is None and "gene" in cl:
            gene_col = c

    # fallback primele 2 coloane
    if drug_col is None or gene_col is None:
        drug_col = df.columns[0]
        gene_col = df.columns[1]

    df = df[[drug_col, gene_col]].dropna()
    df.columns = ["drug", "gene"]

    # construiesc graful bipartit
    G = nx.Graph()
    for drug, gene in df.itertuples(index=False):
        drug = str(drug).strip()
        gene = str(gene).strip()
        if not drug or not gene:
            continue
        G.add_node(drug, bipartite="drug")
        G.add_node(gene, bipartite="gene")
        G.add_edge(drug, gene)

    # ca si drug summary: nr de gene target per drug
    drugs = [n for n, d in G.nodes(data=True) if d.get("bipartite") == "drug"]
    rows = []
    for dname in drugs:
        targets = sorted(G.neighbors(dname))
        rows.append(
            {
                "drug": dname,
                "n_targets": len(targets),
                "targets": ";".join(targets),
            }
        )

    out = pd.DataFrame(rows).sort_values(["n_targets", "drug"], ascending=[False, True])
    out_path = BASE / f"drug_summary_{HANDLE}.csv"
    out.to_csv(out_path, index=False)

    print("[OK] Built bipartite graph")
    print(" - Nodes:", G.number_of_nodes(), "Edges:", G.number_of_edges())
    print("[OK] Wrote:", out_path)


if __name__ == "__main__":
    main()
