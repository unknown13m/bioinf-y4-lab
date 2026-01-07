import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

HANDLE = "unknown13m"
BASE = Path(f"labs/10_integrative/submissions/{HANDLE}")

SNP = BASE / f"snp_matrix_{HANDLE}.csv"
EXPR = BASE / f"expression_matrix_{HANDLE}.csv"
JOINT = BASE / f"multiomics_concat_{HANDLE}.csv"

OUT_SNP = BASE / f"pca_snp_{HANDLE}.png"
OUT_EXPR = BASE / f"pca_expr_{HANDLE}.png"
OUT_JOINT = BASE / f"pca_joint_{HANDLE}.png"


def load_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "Unnamed: 0" in df.columns:
        df = df.set_index("Unnamed: 0")
    else:
        df = df.set_index(df.columns[0])
    df.index = df.index.astype(str).str.strip().str.lower()
    return df


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    keep = []
    for c in out.columns:
        conv = pd.to_numeric(out[c], errors="coerce")
        if conv.notna().any():
            out[c] = conv
            keep.append(c)
    return out[keep]


def ensure_expr_samples_rows(expr: pd.DataFrame, snp_index: pd.Index) -> pd.DataFrame:
    # daca sample-urile sunt pe coloane (index=gene), facem transpose
    cols_lower = pd.Index([str(x).strip().lower() for x in expr.columns])
    idx_lower = pd.Index([str(x).strip().lower() for x in expr.index])
    overlap_cols = snp_index.intersection(cols_lower)
    overlap_idx = snp_index.intersection(idx_lower)

    if len(overlap_cols) > 0 and len(overlap_cols) >= len(overlap_idx):
        expr = expr.copy()
        expr.columns = cols_lower
        expr = expr.T
        expr.index = expr.index.astype(str).str.strip().str.lower()
    else:
        expr = expr.copy()
        expr.index = idx_lower
    return expr


def pca_plot(X: pd.DataFrame, title: str, out_png: Path):
    #PCA 2D
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X.values)

    plt.figure()
    plt.scatter(coords[:, 0], coords[:, 1])

    #etichete cu numele probelor (puține e ok)
    for i, sid in enumerate(X.index):
        plt.annotate(str(sid), (coords[i, 0], coords[i, 1]))

    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print("[OK] saved:", out_png)


def main():
    #SNP-only
    snp = coerce_numeric(load_table(SNP))

    #EXPR-only (rezolvăm orientarea)
    expr_raw = coerce_numeric(load_table(EXPR))
    expr = ensure_expr_samples_rows(expr_raw, snp.index)

    # JOINT (deja e z-score + aliniat)
    joint = load_table(JOINT)

    # aliniem acelasi set de probe pentru SNP si EXPR (sa fie comparabil)
    common = snp.index.intersection(expr.index).intersection(joint.index)
    snp = snp.loc[common]
    expr = expr.loc[common]
    joint = joint.loc[common]

    # standardizare pentru SNP/EXPR aici (JOINT e deja z-score, dar nu strica daca nu re-scalezi)
    snp_scaled = pd.DataFrame(StandardScaler().fit_transform(snp), index=snp.index, columns=snp.columns)
    expr_scaled = pd.DataFrame(StandardScaler().fit_transform(expr), index=expr.index, columns=expr.columns)

    pca_plot(snp_scaled, f"PCA SNP-only ({HANDLE})", OUT_SNP)
    pca_plot(expr_scaled, f"PCA Expression-only ({HANDLE})", OUT_EXPR)
    pca_plot(joint, f"PCA Joint (SNP+Expr) ({HANDLE})", OUT_JOINT)


if __name__ == "__main__":
    main()
