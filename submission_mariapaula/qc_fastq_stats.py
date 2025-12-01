#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
QC simplu pentru FASTQ:
- număr citiri
- lungime medie
- proporție de N
- scor Phred mediu
"""

import argparse
import gzip
from pathlib import Path


def open_maybe_gz(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def qc_fastq(fastq_path: Path):
    total_reads = 0
    total_bases = 0
    n_bases = 0
    phred_sum = 0
    phred_count = 0

    with open_maybe_gz(fastq_path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break  # am ajuns la final
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline().strip()

            if not seq or not qual:
                break  # fișier incomplet

            total_reads += 1
            total_bases += len(seq)
            n_bases += seq.upper().count("N")

            # scoruri Phred (ASCII – 33, format Sanger)
            for ch in qual:
                phred_sum += (ord(ch) - 33)
                phred_count += 1

    if total_reads == 0 or total_bases == 0:
        raise SystemExit("[eroare] FASTQ pare gol sau invalid.")

    avg_len = total_bases / total_reads
    prop_N = n_bases / total_bases
    avg_phred = phred_sum / phred_count if phred_count else 0.0

    return {
        "num_reads": total_reads,
        "avg_length": avg_len,
        "prop_N": prop_N,
        "avg_phred": avg_phred,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True,
                    help="Cale către fișierul FASTQ (poate fi .fastq sau .fastq.gz)")
    ap.add_argument("--out", required=True,
                    help="Fișier text pentru raportul QC")
    args = ap.parse_args()

    fastq_path = Path(args.fastq)
    if not fastq_path.exists():
        raise SystemExit(f"[eroare] Nu găsesc fișierul: {fastq_path}")

    stats = qc_fastq(fastq_path)

    out_path = Path(args.out)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("QC report\n")
        f.write(f"FASTQ: {fastq_path}\n")
        f.write(f"Număr citiri: {stats['num_reads']}\n")
        f.write(f"Lungime medie: {stats['avg_length']:.2f} baze\n")
        f.write(f"Proporție N: {stats['prop_N']:.4f}\n")
        f.write(f"Scor Phred mediu: {stats['avg_phred']:.2f}\n")

    print(f"[OK] Raport QC salvat în {out_path}")


if __name__ == "__main__":
    main()
