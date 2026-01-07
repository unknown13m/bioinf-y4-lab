import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

HANDLE = "unknown13m"
BASE = Path(f"labs/10_integrative/submissions/{HANDLE}")
DATA = BASE / f"multiomics_concat_{HANDLE}.csv"
OUT = BASE / f"cluster_multiomics_{HANDLE}.png"

df = pd.read_csv(DATA, index_col=0)

kmeans = KMeans(n_clusters=2, random_state=42)
labels = kmeans.fit_predict(df)

pca = PCA(n_components=2)
X_pca = pca.fit_transform(df)

plt.figure(figsize=(6,5))
plt.scatter(X_pca[:,0], X_pca[:,1], c=labels, cmap="Set1", s=80)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("KMeans clustering on integrated multi-omics")
plt.tight_layout()
plt.savefig(OUT)
plt.close()

print("[OK] Bonus clustering saved:", OUT)
