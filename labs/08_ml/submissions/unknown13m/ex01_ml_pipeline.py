from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


HANDLE = "unknown13m"
BASE = Path(f"labs/08_ml/submissions/{HANDLE}")


def save_confusion_png(cm, labels, out_png: Path, title: str):
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    ax.imshow(cm)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)

    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title(title)

    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    # primul task 1 aici e partea de load si prepare
  
    in_csv = BASE / f"expression_matrix_{HANDLE}.csv"
    df = pd.read_csv(in_csv)

    if "Condition" not in df.columns:
        raise ValueError("Input CSV must contain a 'Condition' column (labels).")

    y = df["Condition"].astype(str)

    # X sunt gene columns
    drop_cols = [c for c in ["Condition", "SampleID"] if c in df.columns]
    X = df.drop(columns=drop_cols)

    # fac encode labels
    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    # am stratified split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc, test_size=0.30, stratify=y_enc, random_state=42
    )

    # task ul 2 supervised  Random Forest
    rf = RandomForestClassifier(
        n_estimators=300,
        random_state=42,
        class_weight="balanced"
    )
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    # classification report e TXT
    report_txt = classification_report(y_test, y_pred, target_names=le.classes_)
    (BASE / f"classification_report_{HANDLE}.txt").write_text(report_txt)

    # confusion matrix va fi PNG
    cm = confusion_matrix(y_test, y_pred)
    save_confusion_png(
        cm, list(le.classes_),
        BASE / f"confusion_rf_{HANDLE}.png",
        title="Random Forest — Confusion Matrix"
    )

    # feature importance de tip CSV
    fi = pd.DataFrame({
        "gene": X.columns,
        "importance": rf.feature_importances_
    }).sort_values("importance", ascending=False)
    fi.to_csv(BASE / f"feature_importance_{HANDLE}.csv", index=False)

    # al treilea task e logistic regression
    logreg_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=2000))
    ])
    logreg_pipe.fit(X_train, y_train)
    y_pred_lr = logreg_pipe.predict(X_test)

    acc_rf = accuracy_score(y_test, y_pred)
    acc_lr = accuracy_score(y_test, y_pred_lr)

    compare_note = (
        f"Accuracy RF: {acc_rf:.3f}\n"
        f"Accuracy LogReg: {acc_lr:.3f}\n"
        "Note: LogReg is linear + needs scaling; RF is non-linear and handles feature interactions.\n"
    )
    (BASE / f"compare_rf_vs_logreg_{HANDLE}.txt").write_text(compare_note)

    # task ul 4 e unsupervised si face PCA si KMeans
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_pca)

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    ax.scatter(X_pca[:, 0], X_pca[:, 1], c=clusters)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA (2D) + KMeans clusters")
    fig.tight_layout()
    fig.savefig(BASE / f"sup_vs_unsup_scatter_{HANDLE}.png", dpi=200)
    plt.close(fig)

    ctab = pd.crosstab(
        pd.Series(y.values, name="Label"),
        pd.Series(clusters, name="Cluster")
    )
    ctab.to_csv(BASE / f"cluster_crosstab_{HANDLE}.csv")

    # la task ul 5 e semi-supervised (cu pseudo-labeling)
    rng = np.random.RandomState(42)
    y_semi = y_enc.copy()

    mask_unknown = rng.rand(len(y_semi)) < 0.40
    y_semi[mask_unknown] = -1

    rf_semi = RandomForestClassifier(n_estimators=300, random_state=42)
    rf_semi.fit(X[~mask_unknown], y_enc[~mask_unknown])

    pseudo_labels = rf_semi.predict(X[mask_unknown])
    y_full = y_enc.copy()
    y_full[mask_unknown] = pseudo_labels

    rf_pseudo = RandomForestClassifier(n_estimators=300, random_state=42)
    rf_pseudo.fit(X, y_full)

    y_pred_pseudo = rf_pseudo.predict(X_test)
    acc_pseudo = accuracy_score(y_test, y_pred_pseudo)

    semi_note = (
        f"Semi-supervised (pseudo-labeling) with 40% unknown labels\n"
        f"Accuracy RF supervised: {acc_rf:.3f}\n"
        f"Accuracy RF pseudo-labeled: {acc_pseudo:.3f}\n"
        "If pseudo-labels are noisy, performance can drop. If structure is strong, it can improve.\n"
        "Pseudo-labeling helps when labeled data is scarce.\n"
    )
    (BASE / f"semi_supervised_{HANDLE}.txt").write_text(semi_note)

    print("[OK] Done. Outputs written to:", BASE)


if __name__ == "__main__":
    main()
