from __future__ import annotations

from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# aici fac prin NetworkX Louvain
from networkx.algorithms.community import louvain_communities


def load_expression(path: Path) -> pd.DataFrame:
    """
    Încearcă să încarce un tabel de expresie genică.
    Acceptă CSV/TSV.
    Returnează DataFrame cu gene pe rânduri și sample-uri pe coloane.
    """
    # detectez delimiter ul
    if path.suffix.lower() in [".tsv", ".txt"]:
        df = pd.read_csv(path, sep="\t")
    else:
        # aici e fallback csv ul
        df = pd.read_csv(path)

    # dacă e vreo coloană "gene" / "Gene" etc, atunci  o facem index
    gene_col = None
    for c in df.columns:
        if str(c).strip().lower() in ["gene", "genes", "symbol", "gene_symbol", "id"]:
            gene_col = c
            break

    if gene_col is not None:
        df = df.set_index(gene_col)

    # dacă indexul e numeric (0..n) și prima coloană pare gene, dau set index pe prima coloană
    if df.index.dtype.kind in "iu" and df.shape[1] > 1:
        # dacă prima coloană e string uri și restul numeric
        if df.iloc[:, 0].dtype == object:
            df = df.set_index(df.columns[0])

    # păstrez doar numeric
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all").dropna(axis=0, how="all")

    # dacă genele sunt pe coloane și sample urile pe rânduri, fac transpunere (gen heuristic)
    # (mai multe coloane decât rânduri sau indexul pare sample id)
    if df.shape[0] < df.shape[1] and df.shape[0] < 50:
        # posibil rândurile sunt sample uri si coloanele gene, atunci iara transpunem
        df = df.T

    # curat  numele de  gene, eliberez din memorie
    df.index = df.index.astype(str)

    return df


def preprocess(df: pd.DataFrame, min_var_quantile: float = 0.50, max_genes: int = 300) -> pd.DataFrame:
    """
    - log2(x+1)
    - filtrează gene cu varianță scăzută (păstrează peste quantile)
    - limitează nr de gene pentru rețea (ca să fie desenabilă)
    """
    df = df.clip(lower=0)
    df = np.log2(df + 1)

    variances = df.var(axis=1)
    thr = variances.quantile(min_var_quantile)
    df_f = df.loc[variances >= thr].copy()

    # dacă rămân prea multe gene, păstrez topul după varianță
    if df_f.shape[0] > max_genes:
        top = variances.loc[df_f.index].sort_values(ascending=False).head(max_genes).index
        df_f = df_f.loc[top].copy()

    return df_f


def build_network(df: pd.DataFrame, method: str = "pearson", corr_threshold: float = 0.70) -> tuple[nx.Graph, pd.DataFrame]:
    """
    Construiește matrice de corelație gene-gene și un graf (adjacency) cu prag.
    """
    if method == "spearman":
        corr = df.T.corr(method="spearman")
    else:
        corr = df.T.corr(method="pearson")

    # aici daca am zero pe diagonală
    np.fill_diagonal(corr.values, 0.0)

    # adjacency gen daca am  muchii unde corr >= pragul
    mask = corr.abs() >= corr_threshold
    genes = corr.index.tolist()

    G = nx.Graph()
    G.add_nodes_from(genes)

    # adaug muchiile (i<j)
    for i, g1 in enumerate(genes):
        row = mask.iloc[i].values
        for j in range(i + 1, len(genes)):
            if row[j]:
                g2 = genes[j]
                w = float(corr.iloc[i, j])
                G.add_edge(g1, g2, weight=w)

    return G, corr


def run_louvain(G: nx.Graph, seed: int = 42) -> dict[str, int]:
    """
    Detectare module cu Louvain. Returnează dict gene -> module_id (0..k-1).
    """
    if G.number_of_edges() == 0:
        # totul în modul 0
        return {n: 0 for n in G.nodes()}

    comms = louvain_communities(G, seed=seed, weight="weight", resolution=1.0)
    module_map = {}
    for mid, nodes in enumerate(comms):
        for n in nodes:
            module_map[n] = mid
    # noduri izolate care n-au apărut (rar)
    for n in G.nodes():
        module_map.setdefault(n, -1)
    return module_map


def hub_genes(G: nx.Graph, module_map: dict[str, int], top_n: int = 10) -> pd.DataFrame:
    """
    Hub genes = grad (degree) mare (weighted degree = strength) + degree centrality.
    Returnează top N global.
    """
    # strength = suma ponderilor absolute
    strength = {n: 0.0 for n in G.nodes()}
    for u, v, d in G.edges(data=True):
        w = abs(float(d.get("weight", 1.0)))
        strength[u] += w
        strength[v] += w

    deg = dict(G.degree())
    deg_cent = nx.degree_centrality(G) if G.number_of_nodes() > 1 else {n: 0.0 for n in G.nodes()}

    rows = []
    for n in G.nodes():
        rows.append(
            {
                "gene": n,
                "module": int(module_map.get(n, -1)),
                "degree": int(deg.get(n, 0)),
                "strength": float(strength.get(n, 0.0)),
                "degree_centrality": float(deg_cent.get(n, 0.0)),
            }
        )

    df = pd.DataFrame(rows).sort_values(["strength", "degree"], ascending=False).head(top_n)
    return df


def draw_network(G: nx.Graph, module_map: dict[str, int], hubs: set[str], out_png: Path):
    """
    Desenează rețeaua (simplu, pentru un număr moderat de gene).
    Nodi colorați după modul, hubs evidențiate.
    """
    plt.figure(figsize=(12, 10))

    # layout
    pos = nx.spring_layout(G, seed=42, k=None)

    # module ids
    mods = sorted(set(module_map.values()))
    # map module -> color index
    mod_to_idx = {m: i for i, m in enumerate(mods)}

    # culori (colormap discret)
    cmap = plt.get_cmap("tab20")

    node_colors = []
    node_sizes = []
    for n in G.nodes():
        mid = module_map.get(n, -1)
        c = cmap(mod_to_idx.get(mid, 0) % 20)
        node_colors.append(c)
        node_sizes.append(220 if n in hubs else 60)

    # edges (subțiri)
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.6)

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        linewidths=0.8,
        edgecolors="black",
        alpha=0.95,
    )

    # label doar pentru hubs (altfel e haos)
    hub_labels = {h: h for h in hubs if h in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=hub_labels, font_size=8)

    plt.title("TP53 Co-expression Network (colored by Louvain module; hubs labeled)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()


def make_report_pdf(out_pdf: Path, handle: str, chosen_module: int, hubs_df: pd.DataFrame):
    """
    Generează un PDF simplu (max ~2 pagini) cu text + tabel hub genes.
    NOTE: Partea de enrichment (GO/KEGG) trebuie făcută extern (g:Profiler/DAVID).
    Aici includem instrucțiuni și loc pentru rezultate, fără să inventăm termeni.
    """
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    w, h = A4
    c = canvas.Canvas(str(out_pdf), pagesize=A4)

    x = 50
    y = h - 50

    def line(text, dy=14, font="Helvetica", size=11):
        nonlocal y
        c.setFont(font, size)
        c.drawString(x, y, text)
        y -= dy

    line(f"Assignment 6 — Gene Co-Expression Networks (TP53) — {handle}", dy=18, font="Helvetica-Bold", size=14)
    line("")

    line("1) Date si preprocesare", font="Helvetica-Bold")
    line("- log2(x+1) pe expresie")
    line("- filtrare gene cu varianta scazuta (pastrate genele cele mai informative)")
    line("")

    line("2) Retea, module (Louvain)", font="Helvetica-Bold")
    line("- corelatie gene-gene (Pearson/Spearman)")
    line("- prag |corr| pentru adjacency si graf")
    line("- detectie module cu Louvain")
    line("")

    line("3) Hub genes", font="Helvetica-Bold")
    line("Hub genes selectate (top dupa strength/degree):")
    line("gene | module | degree | strength", dy=16, font="Helvetica-Bold", size=10)

    c.setFont("Helvetica", 9)
    for _, r in hubs_df.iterrows():
        if y < 90:
            c.showPage()
            y = h - 50
        c.drawString(x, y, f"{r['gene']} | {int(r['module'])} | {int(r['degree'])} | {r['strength']:.3f}")
        y -= 12

    if y < 140:
        c.showPage()
        y = h - 50

    line("", dy=10)
    line("4) Interpretare biologica si diseaseome", font="Helvetica-Bold")
    line(f"- Modul ales pentru analiza: module {chosen_module}")
    line("- Enrichment (GO/KEGG): rulati g:Profiler / DAVID pe genele din modul ales.")
    line("  In raportul final (manual), lipiti top 5 termeni GO/KEGG + p-values.")
    line("")
    line("Reflectie (diseaseome):", font="Helvetica-Bold")
    line("Acest modul ar putea fi legat de boala daca genele sunt co-activate in acelasi context")
    line("biologic (ex: raspuns la stres, reparare ADN, apoptoza). In diseaseome, hub-urile pot")
    line("actiona ca puncte de convergenta intre fenotipuri si cai moleculare, sugerand subtipuri")
    line("sau potentiale tinte terapeutice.")

    c.showPage()
    c.save()


if __name__ == "__main__":
    HANDLE = "unknown13m"
    base = Path(f"labs/06_wgcna/submissions/{HANDLE}")
    base.mkdir(parents=True, exist_ok=True)

    in_path = base / "tp53.csv"

    out_modules = base / f"modules_tp53_{HANDLE}.csv"
    out_hubs = base / f"hubs_tp53_{HANDLE}.csv"
    out_png = base / f"network_tp53_{HANDLE}.png"
    out_pdf = base / f"report_{HANDLE}.pdf"
    out_zip = base / f"lab6_{HANDLE}.zip"

    # 1) load + preprocess
    df = load_expression(in_path)
    df_p = preprocess(df, min_var_quantile=0.50, max_genes=300)

    # 2) network + modules
    G, corr = build_network(df_p, method="pearson", corr_threshold=0.30)
    module_map = run_louvain(G)

    # export modules (gene -> module)
    pd.DataFrame({"gene": list(module_map.keys()), "module": list(module_map.values())}).to_csv(out_modules, index=False)

    # 3) hubs
    hubs_df = hub_genes(G, module_map, top_n=10)
    hubs_df.to_csv(out_hubs, index=False)

    # pick a module for enrichment: the module of the top hub gene (simple choice)
    chosen_module = int(hubs_df.iloc[0]["module"]) if len(hubs_df) else 0
    hubs_set = set(hubs_df["gene"].tolist())

    # 4) network plot
    draw_network(G, module_map, hubs=hubs_set, out_png=out_png)

    # 5) report pdf (simple)
    make_report_pdf(out_pdf, HANDLE, chosen_module, hubs_df)

    # 6) zip deliverables
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(out_modules, arcname=out_modules.name)
        z.write(out_png, arcname=out_png.name)
        z.write(out_hubs, arcname=out_hubs.name)
        z.write(out_pdf, arcname=out_pdf.name)
        z.write(Path(__file__), arcname=Path(__file__).name)

    print("[OK] Done.")
    print("Input:", in_path)
    print("Edges:", G.number_of_edges(), "Nodes:", G.number_of_nodes())
    print("Saved:", out_modules.name, out_hubs.name, out_png.name, out_pdf.name, out_zip.name)
    print("Chosen module for enrichment:", chosen_module)
