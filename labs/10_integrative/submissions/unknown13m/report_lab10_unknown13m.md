# Lab 10 – Multi-Omics Integration (SNPs + Expression)

##Introducere
Am integrat două straturi omice (SNP-uri și expresie genică)
pentru a identifica relatii cross-omics relevante biologic și pentru a compara analiza single-omics cu analiza integrativă.

##PCA: Single-Omics vs Joint
Am aplicat PCA separat pe:
- matricea SNP,
- matricea de expresie,
- matricea integrată (SNP + expresie).

PCA-ul joint oferă o separare diferită față de analizele single-omics,
unde se vede o variație combinată care nu este vizibilă într-un singur strat omic.

Figuri generate:
- `pca_snp_unknown13m.png`
- `pca_expr_unknown13m.png`
- `pca_joint_unknown13m.png`

##Corelații SNP–Genă
Am calculat corelații Pearson între fiecare SNP și fiecare genă,
păstrând doar perechile cu |r| > 0.5.

Rezultatele sunt salvate în:
- `snp_gene_pairs_unknown13m.csv`

Aceste perechi pot indica potențiale relații funcționale între variația genetică
și nivelul de expresie genică.

##Interpretare biologică
Perechile SNP–genă cu corelație ridicată pot reflecta:
-SNP-uri regulatorii,
-efecte cis sau trans asupra expresiei,
-mecanisme relevante pentru oncologie sau farmacogenomică.

Într-un dataset real, aceste rezultate ar necesita validare experimentală.

##Limitări
- număr mic de probe,
- date de tip demo / sintetice,
- lipsa informațiilor clinice,
- corelația nu implică cauzalitate.

##Concluzie
Integrarea multi-omics oferă informații suplimentare față de analizele
single-omics și permite identificarea relațiilor SNP–genă relevante biologic.
