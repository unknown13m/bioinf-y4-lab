Lab 10 – Multi-Omics Integration (SNPs + Expression)

(Am pastrat structura raportului de la celelalte lab-uri, cum am facut gen la lab 9 cu intro etc)

Introducere
Am integrat doua straturi omice (SNP-uri si expresie genica)
pentru a identifica relatii cross-omics relevante biologic si pentru a compara analiza single-omics cu analiza integrativa.

PCA: Single-Omics vs Joint
Am aplicat PCA separat pe:
- matricea SNP,
- matricea de expresie,
- matricea integrata (SNP + expresie).

PCA-ul joint ofera o separare diferita fata de analizele single-omics,
unde se vede o variatie combinata care nu este vizibila intr-un singur strat omic.

Figuri generate:
-`pca_snp_unknown13m.png`
-`pca_expr_unknown13m.png`
-`pca_joint_unknown13m.png`

Corelatii SNP–Gena
Am calculat corelatii Pearson intre fiecare SNP si fiecare gena,
pastrand doar perechile cu |r| > 0.5.

Rezultatele sunt salvate în `snp_gene_pairs_unknown13m.csv`

Perechile acestea pot indica potentiale relatii functionale intre variatia genetica
si nivelul de expresie genica.

Interpretare biologica
Perechile SNP–gena cu corelatie ridicata pot reflecta:
-SNP-uri regulatorii,
-efecte cis sau trans asupra expresiei,
-mecanisme relevante pentru oncologie sau farmacogenomica.

Intr-un dataset real, aceste rezultate ar necesita validare experimentala.

Limitari
-numar mic de probe,
-date de tip demo / sintetice,
-lipsa informatiilor clinice,
-corelatia nu implica cauzalitate.

Concluzie
Integrarea multi-omics ofera informatii suplimentare fata de analizele
single-omics si permite identificarea relatiilor SNP–gena relevante biologic.
