import pandas as pd
from pathlib import Path
from sklearn.preprocessing import StandardScaler

HANDLE = "unknown13m"
BASE = Path(f"labs/10_integrative/submissions/{HANDLE}")

SNP = BASE / f"snp_matrix_{HANDLE}.csv"
EXPR = BASE / f"expression_matrix_{HANDLE}.csv"
OUT = BASE / f"multiomics_concat_{HANDLE}.csv"


def load_table(path: Path) -> pd.DataFrame:
    """
    aici incarcam CSV-ul si punem prima coloană ca index (de obicei 'Unnamed: 0').
    Nu încearcă să ghicească SampleID din continut (pt că la expression
    prima coloană poate fi GENE, nu sample)
    """
    df = pd.read_csv(path)

    # daca e vreo coloană "Unnamed: 0" (foarte des), o folosim ca index
    if "Unnamed: 0" in df.columns:
        df = df.set_index("Unnamed: 0")
    else:
        # altfel, prima coloană devine index
        df = df.set_index(df.columns[0])

    # sterg gen curatam index-ul
    df.index = df.index.astype(str).str.strip()

    return df


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """
    aici convertim coloanele la numeric unde se poate si daca o coloana nu poate fi convertită, ii o dau drop
    """
    out = df.copy()
    keep_cols = []
    for c in out.columns:
        # încercăm conversie
        converted = pd.to_numeric(out[c], errors="coerce")
        # păstrăm coloana dacă are cel puțin un numeric valid
        if converted.notna().any():
            out[c] = converted
            keep_cols.append(c)
    out = out[keep_cols]
    return out


# 1) load
snp_raw = load_table(SNP)
expr_raw = load_table(EXPR)

# 2) numeric cleanup
snp = coerce_numeric(snp_raw)
expr = coerce_numeric(expr_raw)

# 3) normalizam ID-urile (lowercase) PT SAMPLE-URI
# SNP: sample-urile sunt pe rânduri (index)
snp.index = snp.index.astype(str).str.strip().str.lower()

# expression: poate fi:
#  fie sample-urile pe coloane (index=gene)  si atunci trb transpose
#  sau sample-urile pe randuri (index=sample)si  nu trebuie transpose
# dupa care detectam dupa intersectie daca e overlap intre snp.index si expr.columns, atunci expr trebuie transpose
expr_cols_lower = pd.Index([str(x).strip().lower() for x in expr.columns])
expr_index_lower = pd.Index([str(x).strip().lower() for x in expr.index])

overlap_cols = snp.index.intersection(expr_cols_lower)
overlap_idx = snp.index.intersection(expr_index_lower)

if len(overlap_cols) > 0 and len(overlap_cols) >= len(overlap_idx):
    # cazul clasic daca avem gene pe rânduri, sample-uri pe coloane
    expr = expr.copy()
    expr.columns = expr_cols_lower
    expr = expr.T
    expr.index = expr.index.astype(str).str.strip().str.lower()
else:
    # si cazul in care expresia e deja sample-rows
    expr = expr.copy()
    expr.index = expr_index_lower

# 4)intersectie sample-uri comune
common = snp.index.intersection(expr.index)
print("[INFO] SNP samples:", len(snp), "EXPR samples:", len(expr), "COMMON:", len(common))

if len(common) == 0:
    # debug util
    print("[DEBUG] SNP index head:", list(snp.index[:5]))
    print("[DEBUG] EXPR index head:", list(expr.index[:5]))
    print("[DEBUG] EXPR columns head:", list(expr.columns[:5]))
    raise ValueError("Nu exista sample-uri comune intre SNP si Expression (dupa detect + transpose).")

# 5)aliniem pe sample-uri comune (fix in aceeasi ordine)
common = pd.Index(common)
snp = snp.loc[common]
expr = expr.loc[common]

# 6)Z-score separat pe fiecare omic
snp_z = pd.DataFrame(
    StandardScaler().fit_transform(snp),
    index=snp.index,
    columns=[f"SNP_{c}" for c in snp.columns],
)

expr_z = pd.DataFrame(
    StandardScaler().fit_transform(expr),
    index=expr.index,
    columns=[f"EXP_{c}" for c in expr.columns],
)

# 7) si concatenarea finala
multi = pd.concat([snp_z, expr_z], axis=1)
multi.to_csv(OUT)

print("[OK] wrote:", OUT)
print("[OK] shape:", multi.shape)
