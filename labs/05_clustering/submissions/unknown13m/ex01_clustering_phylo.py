from __future__ import annotations

from pathlib import Path
from collections import Counter
import csv

import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO, Phylo

from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster


def kmer_vector(seq: str, k: int, vocab: list[str]) -> np.ndarray:
    seq = seq.upper().replace("-", "")
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1))
    vec = np.array([counts.get(km, 0) for km in vocab], dtype=float)
    s = vec.sum()
    return vec / s if s > 0 else vec


def build_kmer_matrix(fasta_path: Path, k: int = 3) -> tuple[list[str], np.ndarray]:
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) < 3:
        raise ValueError("Ai nevoie de minim 3 secvențe pentru clustering.")
    names = [r.id for r in records]
    seqs = [str(r.seq) for r in records]

    vocab_set = set()
    for s in seqs:
        s = s.upper().replace("-", "")
        for i in range(len(s) - k + 1):
            vocab_set.add(s[i : i + k])
    vocab = sorted(vocab_set)

    X = np.vstack([kmer_vector(s, k, vocab) for s in seqs])
    return names, X


# aici am facut clustering
def run_kmeans_best_k(X: np.ndarray, k_values=(2, 3, 4, 5), random_state=42):
    best_k = None
    best_score = -1.0
    best_labels = None
    scores = {}

    for k in k_values:
        km = KMeans(n_clusters=k, n_init=10, random_state=random_state)
        labels = km.fit_predict(X)
        if len(set(labels)) < 2:
            continue
        score = silhouette_score(X, labels)
        scores[k] = score
        if score > best_score:
            best_score = score
            best_k = k
            best_labels = labels

    return best_k, best_score, best_labels, scores


def run_hierarchical(X: np.ndarray, method="average", n_clusters=2):
    Z = linkage(X, method=method, metric="euclidean")
    labels = fcluster(Z, t=n_clusters, criterion="maxclust") - 1
    return Z, labels


def run_dbscan(X: np.ndarray, eps=0.5, min_samples=2):
    db = DBSCAN(eps=eps, min_samples=min_samples)
    labels = db.fit_predict(X)  # -1 = noise
    return labels


# un phylo split: 2 clade uri principa;e ale root
def tree_root_split_labels(tree_path: Path, names: list[str]) -> np.ndarray:
    tree = Phylo.read(str(tree_path), "newick")
    root = tree.root

    if len(root.clades) < 2:
        return np.zeros(len(names), dtype=int)

    left = root.clades[0].get_terminals()
    right = root.clades[1].get_terminals()

    left_set = {t.name for t in left}
    right_set = {t.name for t in right}

    labels = []
    for n in names:
        if n in left_set:
            labels.append(0)
        elif n in right_set:
            labels.append(1)
        else:
            labels.append(-1)
    return np.array(labels, dtype=int)


# plot uri
def save_dendrogram(Z, names, out_png: Path):
    plt.figure(figsize=(10, 5))
    dendrogram(Z, labels=names, leaf_rotation=90)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def save_pca_scatter(X, labels, out_png: Path, title: str):
    pca = PCA(n_components=2, random_state=42)
    X2 = pca.fit_transform(X)

    plt.figure(figsize=(7, 5))
    uniq = sorted(set(labels))
    for u in uniq:
        idx = labels == u
        plt.scatter(X2[idx, 0], X2[idx, 1], label=str(u))
    plt.title(title)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def save_silhouette_plot(scores: dict[int, float], out_png: Path):
    ks = sorted(scores.keys())
    vals = [scores[k] for k in ks]
    plt.figure(figsize=(6, 4))
    plt.plot(ks, vals, marker="o")
    plt.title("Silhouette score for KMeans")
    plt.xlabel("K")
    plt.ylabel("Silhouette")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def save_annotated_tree(tree_path: Path, labels: np.ndarray, names: list[str], out_png: Path, prefix: str):
    tree = Phylo.read(str(tree_path), "newick")

    label_map = {n: labels[i] for i, n in enumerate(names)}
    for t in tree.get_terminals():
        if t.name in label_map:
            t.name = f"{t.name}|{prefix}={label_map[t.name]}"

    plt.figure(figsize=(10, 12))
    ax = plt.gca()
    Phylo.draw(tree, do_show=False, axes=ax)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

if __name__ == "__main__":
    base = Path("labs/05_clustering/submissions/unknown13m")
    fasta = base / "input.fasta"
    tree_path = base / "tree_unknown13m.nwk"

    out_labels_csv = base / "cluster_labels.csv"
    out_dendro = base / "dendrogram.png"
    out_pca_kmeans = base / "pca_kmeans.png"
    out_pca_dbscan = base / "pca_dbscan.png"
    out_sil = base / "silhouette_kmeans.png"
    out_tree_kmeans = base / "tree_annotated_kmeans.png"

    # am pus vectorizare
    names, X = build_kmer_matrix(fasta, k=3)

    # apoi KMeans (cel mai bun K cu silhouette)
    best_k, best_sil, kmeans_labels, sil_scores = run_kmeans_best_k(X, k_values=(2, 3, 4, 5))
    if best_k is None:
        raise RuntimeError("KMeans failed to produce valid silhouette scores.")

    save_silhouette_plot(sil_scores, out_sil)
    save_pca_scatter(X, kmeans_labels, out_pca_kmeans, f"PCA (KMeans, best K={best_k}, silhouette={best_sil:.3f})")

    # aici facui Hierarchical clustering (average linkage) + dendrogram
    Z, hier_labels = run_hierarchical(X, method="average", n_clusters=best_k)
    save_dendrogram(Z, names, out_dendro)

    # DBSCAN 
    dbscan_labels = run_dbscan(X, eps=0.5, min_samples=2)
    save_pca_scatter(X, dbscan_labels, out_pca_dbscan, "PCA (DBSCAN, -1=noise)")

    # am comparat cu un phylogenetic root split (2 clade uri)
    phylo_labels = tree_root_split_labels(tree_path, names)

    valid = phylo_labels != -1
    ari_kmeans = adjusted_rand_score(phylo_labels[valid], kmeans_labels[valid])
    ari_hier = adjusted_rand_score(phylo_labels[valid], hier_labels[valid])

    save_annotated_tree(tree_path, kmeans_labels, names, out_tree_kmeans, prefix=f"KMeansK{best_k}")

    with out_labels_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sequence", "phylo_root_clade", f"kmeans_best_k_{best_k}", "hierarchical", "dbscan"])
        for i, n in enumerate(names):
            w.writerow([n, int(phylo_labels[i]), int(kmeans_labels[i]), int(hier_labels[i]), int(dbscan_labels[i])])

    print("[OK] Outputs generated in:", base)
    print(" - dendrogram:", out_dendro)
    print(" - PCA kmeans:", out_pca_kmeans)
    print(" - PCA dbscan:", out_pca_dbscan)
    print(" - silhouette:", out_sil)
    print(" - tree annotated:", out_tree_kmeans)
    print(" - labels csv:", out_labels_csv)
    print("\nComparison vs phylogenetic root split (Adjusted Rand Index):")
    print(f" - ARI KMeans(best K={best_k}): {ari_kmeans:.3f}")
    print(f" - ARI Hierarchical(n_clusters={best_k}): {ari_hier:.3f}")
