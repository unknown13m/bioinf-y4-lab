from pathlib import Path

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

if __name__ == "__main__":
    
    fasta = Path("labs/04_phylogenetics/submissions/unknown13m/input.fasta")

    # aici citesc alinierea (Biopython vede fasta ca alignment daca is lungimi similare)
    alignment = AlignIO.read(fasta, "fasta")

    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    output = Path("labs/04_phylogenetics/submissions/unknown13m/tree_unknown13m.nwk")
    output.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, output, "newick")

    print("\n[OK] Neighbor-Joining tree (ASCII):\n")
    Phylo.draw_ascii(tree)

    print(f"\n[OK] Saved Newick to: {output}\n")

    try:
        import matplotlib.pyplot as plt

        png = Path("labs/04_phylogenetics/submissions/unknown13m/tree_unknown13m.png")
        fig = plt.figure(figsize=(10, 12))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, axes=ax)
        fig.tight_layout()
        fig.savefig(png, dpi=200)
        plt.close(fig)
        print(f"[OK] Saved PNG to: {png}\n")
    except Exception as e:
        print(f"[INFO] PNG not generated (matplotlib missing or error): {e}")
