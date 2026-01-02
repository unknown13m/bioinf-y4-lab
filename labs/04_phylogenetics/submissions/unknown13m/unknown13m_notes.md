Lab4 — Filogenetică (distanțe + arbore NJ)

Ce secvențe FASTA am folosit:
-un fișier multi-FASTA cu 10 secvențe ADN:
`labs/04_phylogenetics/submissions/unknown13m/input.fasta`.

Setul meu are două grupuri :
-clusterA: seq01_clusterA – seq05_clusterA;
-clusterB: seq06_clusterB – seq10_clusterB.

Matricea de distanțe: am calculat matricea de distanțe și am salvat-o în `labs/04_phylogenetics/submissions/unknown13m/distance_matrix.csv`.

Observații:
-distanțele între secvențele din același cluster sunt mici (secvențe foarte similare);
-distanțele dintre clusterA și clusterB sunt mai mari (secvențe mai diferite).

La arbore Neighbor-Joining (Biopython):
-aici am construit un arbore Neighbor-Joining cu Biopython (DistanceCalculator + DistanceTreeConstructor);
Si dupa am salvat arborele în format Newick:
`labs/04_phylogenetics/submissions/unknown13m/tree_unknown13m.nwk`

in plus, am generat și o vizualizare grafică:
`labs/04_phylogenetics/submissions/unknown13m/tree_unknown13m.png`.

interpretare clustere:
-seq01_clusterA – seq05_clusterA se grupează împreună (clusterA);
-seq06_clusterB – seq10_clusterB se grupează împreună (clusterB);
-cele două clustere sunt separate si de asta am două grupuri distincte în funcție de similaritatea secvențelor.

MSA (Clustal Omega) și regiuni conservate:
-dacă introduc multi-FASTA-ul în Clustal Omega, în aliniere se observă:
1.multe poziții conservate în interiorul clusterA (coloane aproape identice);
2.multe poziții conservate în interiorul lui clusterB;
3.mai puține poziții conservate global când compar clusterA cu clusterB (diferențe pe multe poziții).

Concluzie: rezultatul MSA este consistent cu arborele NJ — secvențele cu regiuni conservate similare apar în același cluster.


Astfel, arborele filogenetic oferă informație suplimentară față de o matrice de distanțe deoarece:
-oferă o reprezentare ierarhică a înrudirii (clustere/subclustere);
-face vizibil rapid cine este mai apropiat de cine;
-este mai ușor de interpretat decât doar numerele dintr-un tabel de distanțe.
