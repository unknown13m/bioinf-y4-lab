import pandas as pd
import numpy as np
from pathlib import Path

HANDLE = "unknown13m"
BASE = Path(f"labs/10_integrative/submissions/{HANDLE}")

SNP = BASE / f"snp_matrix_{HANDLE}.csv"
EXPR = BASE / f"expression_matrix_{HANDLE}.csv"
OUT = BASE / f"snp_gene_pairs_{HANDLE}.csv"


def load_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "Unnamed: 0" in df.columns:
        df = df.set_index("Unnamed: 0")
    else:
        df = df.set_index(df.columns[0])
    df.index = df.index.astype(str).str.strip().str.lower()
    return df


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    out = {}
    for c in df.columns:
        v = pd.to_numeric(df[c], errors="coerce")
        if v.notna().sum() > 1:
            out[c] = v
    return pd.DataFrame(out, index=df.index)


def ensure_expr_samples_rows(expr: pd.DataFrame, snp_index: pd.Index) -> pd.DataFrame:
    cols_lower = pd.Index([str(x).strip().lower() for x in expr.columns])
    idx_lower = pd.Index([str(x).strip().lower() for x in expr.index])

    if len(snp_index.intersection(cols_lower)) >= len(snp_index.intersection(idx_lower)):
        expr = expr.copy()
        expr.columns = cols_lower
        expr = expr.T
        expr.index = expr.index.astype(str).str.strip().str.lower()
    else:
        expr.index = idx_lower
    return expr


def main():
    snp = coerce_numeric(load_table(SNP))
    expr_raw = coerce_numeric(load_table(EXPR))
    expr = ensure_expr_samples_rows(expr_raw, snp.index)

    common = snp.index.intersection(expr.index)
    snp = snp.loc[common]
    expr = expr.loc[common]

    rows = []
    for snp_col in snp.columns:
        x = snp[snp_col].values
        for gene_col in expr.columns:
            y = expr[gene_col].values
            if np.std(x) == 0 or np.std(y) == 0:
                continue
            r = np.corrcoef(x, y)[0, 1]
            if abs(r) > 0.5:
                rows.append({
                    "snp": snp_col,
                    "gene": gene_col,
                    "pearson_r": round(r, 3)
                })

    out = pd.DataFrame(rows).sort_values("pearson_r", key=lambda s: s.abs(), ascending=False)
    out.to_csv(OUT, index=False)
    print("[OK] wrote:", OUT, "pairs:", len(out))


if __name__ == "__main__":
    main()
