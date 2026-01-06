import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path

HANDLE = "unknown13m"
BASE = Path(f"labs/09_repurposing/submissions/{HANDLE}")
DATA = Path("data/sample/drug_gene_tp53_subset.csv")

OUT_PNG = BASE / f"network_drug_gene_{HANDLE}.png"

def main():
    BASE.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    # expected columns: drug,gene,... (we only need drug + gene)
    if not {"drug", "gene"}.issubset(df.columns):
        raise ValueError(f"Expected columns drug,gene. Got: {df.columns.tolist()}")

    # Build bipartite graph
    G = nx.Graph()
    drugs = sorted(df["drug"].dropna().unique().tolist())
    genes = sorted(df["gene"].dropna().unique().tolist())

    for d in drugs:
        G.add_node(d, bipartite="drug")
    for g in genes:
        G.add_node(g, bipartite="gene")

    for _, row in df.iterrows():
        d = row["drug"]
        g = row["gene"]
        if pd.isna(d) or pd.isna(g):
            continue
        G.add_edge(d, g)

    # Node sizes = number of gene targets (for drugs) / degree (for genes)
    drug_targets = df.groupby("drug")["gene"].nunique().to_dict()
    sizes = []
    colors = []
    for n in G.nodes():
        if n in drug_targets:
            # drugs
            colors.append("tab:blue")
            sizes.append(600 + 600 * drug_targets.get(n, 1))  # scale
        else:
            # genes
            colors.append("tab:red")
            sizes.append(500 + 200 * G.degree(n))  # scale

    # Layout
    pos = nx.spring_layout(G, seed=42, k=1.2)

    plt.figure(figsize=(12, 9))
    nx.draw_networkx_edges(G, pos, alpha=0.35, width=1.2)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=sizes, linewidths=0.8)
    nx.draw_networkx_labels(G, pos, font_size=8)

    plt.title(f"Drug–Gene network ({HANDLE})\nblue=drugs, red=genes; node size ~ #targets/degree")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=200)
    plt.close()

    print("[OK] Saved:", OUT_PNG)
    print(f"[OK] Nodes: {G.number_of_nodes()} Edges: {G.number_of_edges()}")

if __name__ == "__main__":
    main()
