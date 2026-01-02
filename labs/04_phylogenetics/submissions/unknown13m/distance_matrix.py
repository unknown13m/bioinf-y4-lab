from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import List, Tuple

from Bio import SeqIO


def p_distance(a: str, b: str) -> float:
    """    gen aici ignor gap-uri si ce nu i litera
    """
    a = a.upper()
    b = b.upper()
    n = 0
    diff = 0
    L = min(len(a), len(b))
    for i in range(L):
        ca = a[i]
        cb = b[i]
        if ca == "-" or cb == "-":
            continue
        if not ca.isalpha() or not cb.isalpha():
            continue
        n += 1
        if ca != cb:
            diff += 1
    return (diff / n) if n > 0 else 0.0


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records = [(r.id, str(r.seq)) for r in SeqIO.parse(str(path), "fasta")]
    if len(records) < 10:
        raise ValueError(f"Trebuie >=10 secvențe, dar ai {len(records)} în {path}")
    return records


def build_matrix(records: List[Tuple[str, str]]) -> Tuple[List[str], List[List[float]]]:
    names = [rid for rid, _ in records]
    seqs = [seq for _, seq in records]
    n = len(seqs)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = p_distance(seqs[i], seqs[j])
            mat[i][j] = d
            mat[j][i] = d
    return names, mat


def write_csv(out_csv: Path, names: List[str], mat: List[List[float]]) -> None:
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + names)
        for name, row in zip(names, mat):
            w.writerow([name] + [f"{x:.6f}" for x in row])


def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: python distance_matrix.py <input.fasta> <out.csv>")
        sys.exit(1)

    in_fa = Path(sys.argv[1])
    out_csv = Path(sys.argv[2])

    records = read_fasta(in_fa)
    names, mat = build_matrix(records)
    write_csv(out_csv, names, mat)

    print(f"[OK] Distance matrix saved to: {out_csv}")


if __name__ == "__main__":
    main()
