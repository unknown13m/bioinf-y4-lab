import gzip
import matplotlib.pyplot as plt
from collections import Counter

def read_fastq(path):
    lengths = []
    phreds = []

    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()

            lengths.append(len(seq))
            phreds.extend([ord(q) - 33 for q in qual])

    return lengths, phreds


if __name__ == "__main__":
    fastq = "your_reads.fastq"
    lengths, phreds = read_fastq(fastq)

    plt.figure(figsize=(10,5))

    plt.subplot(1,2,1)
    plt.hist(lengths, bins=20)
    plt.title("Distribuția lungimilor citirilor")
    plt.xlabel("Lungime")
    plt.ylabel("Număr citiri")

    plt.subplot(1,2,2)
    plt.hist(phreds, bins=20)
    plt.title("Distribuția scorurilor Phred")
    plt.xlabel("Scor Phred")
    plt.ylabel("Frecvență")

    plt.tight_layout()
    plt.savefig("qc_plot_mariapaula.png")
    print("[OK] Am salvat graficul în qc_plot_mariapaula.png")
