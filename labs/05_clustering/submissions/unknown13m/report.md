Lab 5 — Clustering și Filogenetică

Pentru dataset-ul și arborele refolosit din Lab4 am refolosit:
-multi-FASTA: `input.fasta`
-arborele filogenetic: `tree_unknown13m.nwk`

Dataset-ul conține 10 secvențe ADN (două grupuri ), iar arborele din Lab4 separă secvențele în două clade principale.

La partea de preprocesare am reprezentat fiecare secvență ca vector de caracteristici folosind frecvențe normalizate de k-mer (k=3), asa ca permite clustering și proiecție PCA în 2D.

La clustering am facut minim 2 metode + evaluarea
aici e KMeans: am testat K ∈ {2,3,4,5} și am ales K optim pe baza Silhouette Score (grafic în `silhouette_kmeans.png`), iar rezultatul e vizualizat în `pca_kmeans.png`.

Apoi la hierarchical (average linkage) am folosit average linkage și am generat dendrograma în `dendrogram.png`, iar nr de clustere folosit este același cu best-K din KMeans pentru comparabilitate.

La DBSCAN am rulat DBSCAN (eps=0.5, min_samples=2), iar rezultatul e în `pca_dbscan.png` (label -1 e noise ul).

A partea de integrare cu filogenetica am comparat clusterele cu arborele asa:
-inatai am separat arborele în 2 clade principale (root split),
-si dupa am comparat etichetele de clustering cu etichetele din arbore.

Apoi am facut concordanța cu Adjusted Rand Index (ARI) ce s a afisat în consola.
Arborele cu etichetele KMeans este in `tree_annotated_kmeans.png`.


-KMeans și hierarchical vor să recupereze separarea în două grupuri observată și în arbore
-DBSCAN poate produce rezultate diferite în funcție de parametrii eps/min_samples (sensibil pe seturi mici)
-diferențele apar deoarece arborele are o structură evolutivă, iar clustering-ul depinde de reprezentarea numerică (k-mer) și de metoda statistică pe care o aleg

1)Cum se aliniază clustering-ul cu arborele?  
În general, clusterele obținute sunt consistente cu separarea în două clade principale.

2)Ce explică diferențele între metode și tree vs clustering?  
Metodele optimizează criterii diferite (densitate, centroid, ierarhie), iar tree araata distanțe/structură evolutivă, reprezentarea (k-mer) influențează clustering-ul.

3)Utilitatea combinării clustering + filogenetică  
Validarea grupărilor, identificarea subtipurilor, familii de gene, posibile diferențe funcționale și ipoteze despre evoluție (diversificare / convergență).
