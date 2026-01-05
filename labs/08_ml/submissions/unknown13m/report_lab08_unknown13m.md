Lab 08: Machine Learning pe date omice


Am folosit un set de date de tip “expression matrix” cu etichete (Condition) și variabile (gene/features).
am facut:
1.clasificare supervised (Random Forest),
2.explorare unsupervised (PCA + KMeans) si
un pseudo-labeling de tipul semi-supervised.

datele mele sunt in:`expression_matrix_unknown13m.csv`  
-y (label): `Condition`
-X (features): Gene001…Gene050 (aici am valori numerice)

si ca si limitare importantă: dataset-ul este mic, deci rezultatele pot varia și nu sunt stabile statistic.

La partea de supervised ML, Random Forest,  am antrenat un RandomForestClassifier (n_estimators=300) pe split stratificat train/test (70/30) si am generat:
-`classification_report_unknown13m.txt`
-`confusion_rf_unknown13m.png`
-`feature_importance_unknown13m.csv`

am vazut ca la confusion matrix arată unde modelul greșește între clase si ca la feature importances imi da un top al genelor cu impact mai mare în decizie.

La logistic regression am rulat și logistic regression cu scaling (StandardScaler) și am comparat cu RF, iar rezultatele sunt in `compare_rf_vs_logreg_unknown13m.txt`.

ca si diferență de concept, la LogReg am un model linear (decizie liniară în spațiul feature-urilor), iar la RF am un model non-linear (captură interacțiuni între feature-uri).

La unsupervised ML si PCA + KMeans am aplicat StandardScaler + PCA (2 componente) pentru vizualizare, apoi KMeans (K=2), dupa care am salvat  scatter PCA colorat după cluster: `sup_vs_unsup_scatter_unknown13m.png` si tabel Label x Cluster: `cluster_crosstab_unknown13m.csv`


-dacă etichetele (Condition) se aliniază cu clusterele, atunci structura nesupravegheată recuperează parțial clasele
-probe “greșite” pot apărea din zgomot, outlieri sau lipsă de separare

Unde am avut semi-supervised (pseudo-labeling) am simulat lipsa etichetelor: ~40% etichete marcate ca “unknown” (-1).
intai am facut RF antrenat doar pe date etichetate, dupa care prezic pseudo-labels pentru cele fără etichetă si la final fac reantrenare pe set complet (etichete reale + pseudo).

Rezumatul e in `semi_supervised_unknown13m.txt`.


-performanța poate crește dacă pseudo-labels sunt corecte și structura e clară, dar poate si scădea dacă pseudo-labels introduc erori (propagare de zgomot).

ca si interpretare biologică (la nivel de lab), cele mai importante “gene” (features) sunt cele cu importanță mare în RF, iar intr-un dataset real, acestea ar putea corespunde unor gene marker asociate cu condiția (Tumor/Normal).
Doar ca aici, deoarece datele sunt de demo și dataset-ul este mic, interpretarea biologică trebuie tratată doar ca exemplu tocmai din cauza datelor putine.

Am avut si limitări de genul nr mic de probe (si asta mi-a cam dat rezultate instabile), un posibil bias din modul de generare a datelor sau din distribuții si imi lipseste variabilitatea biologică (tocmai fiindca am un dataset demo).
