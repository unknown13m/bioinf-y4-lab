from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm

HANDLE = "unknown13m"
BASE = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_PDF = BASE / f"report_repurposing_{HANDLE}.pdf"
NET_IMG = BASE / f"network_drug_gene_{HANDLE}.png"

def main():
    c = canvas.Canvas(str(OUT_PDF), pagesize=A4)
    w, h = A4
    x = 2 * cm
    y = h - 2 * cm
    lh = 14

    c.setFont("Helvetica-Bold", 14)
    c.drawString(x, y, "Lab 09 — Drug Repurposing (Network-Based Approaches)")
    y -= 2 * lh

    c.setFont("Helvetica", 11)
    text = [
        "1. Metodologie",
        "- Am construit o retea bipartita drug–gene folosind dataset-ul drug_gene_tp53_subset.csv.",
        "- Am proiectat reteaua pe layer-ul de medicamente folosind similaritatea Jaccard.",
        "- Am calculat proximitatea medicamentelor fata de genele bolii in reteaua drug–gene.",
        "",
        "2. Cele mai apropiate medicamente",
        "- Medicamentele cu proximitate minima fata de genele bolii sunt prioritizate.",
        "- Exemple observate: Cisplatin, Doxorubicin, Etoposide.",
        "",
        "3. Interpretare biologica",
        "- Medicamentele identificate actioneaza asupra genei TP53 sau a regulatorilor sai (ex: MDM2).",
        "- Acest lucru este consistent cu mecanisme cunoscute de raspuns la DNA damage si cancer.",
        "",
        "4. Limitari",
        "- Dataset-ul este unul de dimensiune redusa si demonstrativ.",
        "- Reteaua este partiala si depinde de adnotarile disponibile.",
        "- Rezultatele trebuie interpretate ca exemplu metodologic, nu concluzie clinica."
    ]

    for line in text:
        if y < 2 * cm:
            c.showPage()
            c.setFont("Helvetica", 11)
            y = h - 2 * cm
        c.drawString(x, y, line)
        y -= lh

    # Pagina cu figura
    c.showPage()
    c.setFont("Helvetica-Bold", 12)
    c.drawString(x, h - 2 * cm, "Vizualizare retea drug–gene")

    if NET_IMG.exists():
        c.drawImage(str(NET_IMG), x, h - 14 * cm, width=16 * cm, preserveAspectRatio=True)

    c.save()
    print("[OK] PDF generated:", OUT_PDF)

if __name__ == "__main__":
    main()
