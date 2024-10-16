# Explications des fonctions

## Fichier `hash.py`

### `filter_smallest(seq, k, s)`
Cette fonction parcourt une séquence de k-mers et sélectionne les `s` plus petits k-mers en utilisant une structure de tas (`heapq`). Les k-mers sont stockés sous forme négative pour simuler un max-heap, ce qui permet de garder les plus petits éléments.

### `hyperMinSketch(seq1, seq2)`
Cette fonction compare deux listes triées de k-mers (séquences) et calcule l'intersection et l'union entre elles. Cela permet de calculer des mesures de similarité comme la similarité de Jaccard en se basant sur les k-mers communs.

### `xorshift(val)`
Un générateur de nombres pseudo-aléatoires basé sur l'algorithme **xorshift**. Cet algorithme manipule les bits d'un entier (`val`) pour générer des valeurs aléatoires utilisées dans le hachage des k-mers.

### `partition_min_hash(kmer_list, num_partitions=1000, num_hashes=1000)`
Cette fonction divise une liste de k-mers en partitions (buckets), puis génère des signatures MinHash pour chaque partition. Chaque partition est hachée plusieurs fois avec différentes fonctions de hachage pour générer des signatures minimales qui sont ensuite utilisées pour comparer les séquences.

---

## Fichier `__main__.py`

### `jaccard_hyperMinSketch_from_files(file1, file2, k, s)`
Cette fonction charge deux fichiers de séquences FASTA, applique un filtrage pour sélectionner les plus petits k-mers, puis utilise la fonction `hyperMinSketch` pour calculer l'intersection et l'union des k-mers entre les deux fichiers. Enfin, elle retourne la similarité de Jaccard en divisant l'intersection par l'union.

### `jaccard_partition_minhash(seq1, seq2, num_partitions, num_hashes)`
Cette fonction calcule la similarité de Jaccard entre deux séquences en utilisant la technique de **Partition MinHash**. Elle génère des signatures MinHash pour chaque partition de k-mers, puis compare ces signatures pour calculer l'intersection et l'union, et retourne la similarité de Jaccard.

Il y a des problèmes cela ne marche pas. J'obtiens 0. Je n'ai pas réussi à télécharger les données pour Homme,souris et singe. J'avais des fichiers vides.

