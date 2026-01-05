import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from pathlib import Path

HANDLE = "unknown13m"
BASE = Path(f"labs/08_ml/submissions/{HANDLE}")

# incep cu load data
df = pd.read_csv(BASE / f"expression_matrix_{HANDLE}.csv")

X = df.drop(columns=["SampleID", "Condition"])
y = df["Condition"]

# aici am PCA inainte de filtering
X_scaled = StandardScaler().fit_transform(X)
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

plt.figure()
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=(y == "Tumor"), alpha=0.7)
plt.title("PCA before variance filtering")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig(BASE / f"pca_before_variance_{HANDLE}.png")
plt.close()

# filtrarea variantei
variances = X.var()
low_var_genes = variances.sort_values().head(10).index
X_filt = X.drop(columns=low_var_genes)

# PCA AFTER filtering
Xf_scaled = StandardScaler().fit_transform(X_filt)
Xf_pca = pca.fit_transform(Xf_scaled)

plt.figure()
plt.scatter(Xf_pca[:, 0], Xf_pca[:, 1], c=(y == "Tumor"), alpha=0.7)
plt.title("PCA after removing low-variance genes")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig(BASE / f"pca_after_variance_{HANDLE}.png")
plt.close()

print("[OK] Bonus PCA variance filtering done.")
print("Removed genes:", list(low_var_genes))
