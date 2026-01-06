import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path

HANDLE = "unknown13m"

BASE = Path(f"labs/09_repurposing/submissions/{HANDLE}")
DATA = Path("data/sample/drug_gene_tp53_subset.csv")
DISEASE_GENES = BASE / f"disease_genes_{HANDLE}.txt"

OUT_RWR = BASE / f"drug_priority_rwr_{HANDLE}.csv"
OUT_COMPARE = BASE / f"compare_proximity_vs_rwr_{HANDLE}.csv"

RESTART = 0.5
MAX_ITERS = 500
TOL = 1e-10


def build_bipartite_graph(df: pd.DataFrame) -> nx.Graph:
    """
    Build drug-gene bipartite graph (undirected).
    Nodes are raw names from CSV.
    """
    G = nx.Graph()
    for _, r in df.iterrows():
        drug = str(r["drug"]).strip()
        gene = str(r["gene"]).strip()
        if not drug or not gene:
            continue
        G.add_node(drug, kind="drug")
        G.add_node(gene, kind="gene")
        G.add_edge(drug, gene)
    return G


def rwr_scores(G: nx.Graph, seed_nodes: list[str], restart: float = 0.5) -> pd.Series:
    """
    Random Walk with Restart (numpy version, no scipy):
      p_{t+1} = (1-restart) * W * p_t + restart * p0
    W is column-stochastic adjacency.
    """
    nodes = list(G.nodes())
    n = len(nodes)
    idx = {node: i for i, node in enumerate(nodes)}

    seeds = [s for s in seed_nodes if s in idx]
    if len(seeds) == 0:
        return pd.Series(0.0, index=nodes)

    # adiacenta  (numpy)
    A = nx.to_numpy_array(G, nodelist=nodes, dtype=float)

    # normalizare column-stochastic 
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    W = A / col_sums  # broadcast imparte fiecare coloana

    # p0
    p0 = np.zeros(n, dtype=float)
    for s in seeds:
        p0[idx[s]] = 1.0
    p0 /= p0.sum()

    p = p0.copy()
    for _ in range(MAX_ITERS):
        p_next = (1.0 - restart) * (W @ p) + restart * p0
        if np.linalg.norm(p_next - p, ord=1) < TOL:
            p = p_next
            break
        p = p_next

    return pd.Series(p, index=nodes)



def main():
    BASE.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    # coloane gen drug,gene,mechanism,source_note (am nevoie doar de drug,gene)
    if not {"drug", "gene"}.issubset(df.columns):
        raise ValueError(f"CSV must contain columns drug,gene. Found: {df.columns.tolist()}")

    G = build_bipartite_graph(df)

    # citeste disease genes seeds
    if not DISEASE_GENES.exists():
        raise FileNotFoundError(f"Missing disease genes file: {DISEASE_GENES}")

    seeds = [x.strip() for x in DISEASE_GENES.read_text(encoding="utf-8").splitlines() if x.strip()]
    # le pastrez doar pe cele din graf
    seeds_in_graph = [s for s in seeds if s in G]

    scores = rwr_scores(G, seeds_in_graph, restart=RESTART)

    # facem scorul pt drugs
    drugs = [n for n, d in G.nodes(data=True) if d.get("kind") == "drug"]
    drug_scores = scores.loc[drugs].sort_values(ascending=False)

    out = pd.DataFrame({
        "drug": drug_scores.index,
        "rwr_score": drug_scores.values
    })
    out.to_csv(OUT_RWR, index=False)

    #daca am deja un output care sa ne fie cata pe aproape de ce trebuie, compar ranking ul
    prox_path = BASE / f"drug_priority_{HANDLE}.csv"
    if prox_path.exists():
        prox = pd.read_csv(prox_path)
        prox = prox.rename(columns={"mean_distance": "proximity_mean_distance"})
        merged = prox.merge(out, on="drug", how="outer")
        merged.to_csv(OUT_COMPARE, index=False)

        print("[OK] Wrote:", OUT_RWR)
        print("[OK] Wrote:", OUT_COMPARE)
    else:
        print("[OK] Wrote:", OUT_RWR)
        print("[WARN] Proximity file not found, skipping compare:", prox_path)


if __name__ == "__main__":
    main()
