Lab 08: Machine Learning pe date omice

1. **INTRODUCERE**
Am folosit un set de date de tip “expression matrix” cu etichete (Condition) si variabile (gene/features).
Am facut:
1.clasificare supervised (Random Forest),
2.explorare unsupervised (PCA + KMeans) si
un pseudo-labeling de tipul semi-supervised.

Datele mele sunt in:`expression_matrix_unknown13m.csv`  
-y (label): `Condition`
-X (features): Gene001…Gene050 (aici am valori numerice)

Iar ca si limitare importanta: dataset-ul este mic, deci rezultatele pot varia si nu sunt stabile statistic.

2. **SUPERVISED ML (Random Forest/Logistic Regression)**
La partea de supervised ML, Random Forest,  am antrenat un RandomForestClassifier (n_estimators=300) pe split stratificat train/test (70/30) si am generat:
-`classification_report_unknown13m.txt`
-`confusion_rf_unknown13m.png`
-`feature_importance_unknown13m.csv`

Am vazut ca la confusion matrix arata unde modelul greseste intre clase si ca la feature importances imi da un top al genelor cu impact mai mare in decizie.

La logistic regression am rulat si logistic regression cu scaling (StandardScaler) si am comparat cu RF, iar rezultatele sunt in `compare_rf_vs_logreg_unknown13m.txt`.

Ca si diferenta de concept, la LogReg am un model linear (decizie liniara in spatiul feature-urilor), iar la RF am un model non-linear (captura interactiuni intre feature-uri).

3. **UNSUPERVISED ML (PCA + KMeans)**
La unsupervised ML si PCA + KMeans am aplicat StandardScaler + PCA (2 componente) pentru vizualizare, apoi KMeans (K=2), dupa care am salvat  scatter PCA colorat dupa cluster: `sup_vs_unsup_scatter_unknown13m.png` si tabel Label x Cluster: `cluster_crosstab_unknown13m.csv`


-daca etichetele (Condition) se aliniaza cu clusterele, atunci structura nesupravegheata recupereaza partial clasele
-probe “gresite” pot aparea din zgomot, outlieri sau lipsa de separare


4. **SEMI-SUPERVISED LEARNING**
Unde am avut semi-supervised (pseudo-labeling) am simulat lipsa etichetelor: ~40% etichete marcate ca “unknown” (-1).
Intai am facut RF antrenat doar pe date etichetate, dupa care prezic pseudo-labels pentru cele fara eticheta si la final fac reantrenare pe set complet (etichete reale + pseudo).

Rezumatul e in `semi_supervised_unknown13m.txt`.


-performanta poate creste daca pseudo-labels sunt corecte si structura e clara, dar poate si scadea daca pseudo-labels introduc erori (propagare de zgomot).


5. **INTERPRETARE BIOLOGICA**
Ca si interpretare biologica (la nivel de lab), cele mai importante “gene” (features) sunt cele cu importanta mare în RF, iar intr-un dataset real, acestea ar putea corespunde unor gene marker asociate cu conditia (Tumor/Normal).
Doar ca aici, deoarece datele sunt de demo si dataset-ul este mic, interpretarea biologica trebuie tratata doar ca exemplu tocmai din cauza datelor putine.


6. **LIMITARI**
Am avut si limitari de genul nr mic de probe (si asta mi-a cam dat rezultate instabile), un posibil bias din modul de generare a datelor sau din distributii si imi lipseste variabilitatea biologica (tocmai fiindca am un dataset demo).


7. **BONUS**
La partea de bonus cu PCA si filtrare dupa varianta, dupa eliminarea a 10 gene cu varianta scazuta, proiectia PCA are o separare putin mai clara intre probe si filtrarea reduce zgomotul care venea din feature-uri neinformative.
Doar ca datorita dimensiunii mici a dataset-ului, efectul este moderat.
