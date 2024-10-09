# Rapport de rendu 

AZMOUDEH Anastasia

L'objectif de ce TP est de comparer un ensemble de 5 bactéries pour retrouver les familles présentes. L'objectif final est de produire une matrice des indices de Jaccard entre tous les échantillons


## Méthodes implémentées

### Fichier `kmers.py`

#### `kmer2str(val, k)`
Transforme un entier représentant un k-mer en sa représentation sous forme de chaîne de caractères. Cette fonction est utile pour convertir les k-mers encodés en chaînes lisibles.

#### `encode_nucl(lettre)`
Encode un nucléotide (A, C, T, G) sous forme d'entier (A=0, C=1, T=2, G=3). Cette fonction permet de transformer une séquence en une représentation numérique.

#### `encode_kmer(seq, k)`
Encode une séquence k-mer en entier. Cela permet d'optimiser le stockage et le traitement des k-mers.

#### `enumerate_kmers(seq, k)`
Génère tous les k-mers d'une séquence donnée en utilisant leur forme encodée en entier. Cette fonction est cruciale pour l'extraction des k-mers des séquences ADN.

#### `reverse_complement(kmer)`
Renvoie le complément inverse d'un k-mer. Cette méthode est utile pour gérer les séquences ADN en double brin où un k-mer et son complément inverse peuvent être considérés comme équivalents.

#### `stream_kmers(text, k)`
Renvoie les k-mers canoniques (le plus petit lexicographiquement entre un k-mer et son complément inverse). Cela assure que chaque k-mer est traité de manière cohérente.

### Fichier `__main__.py`

#### `jaccard(fileA, fileB, k)`
Calcule la similarité de Jaccard entre les k-mers de deux séquences ADN. Les k-mers uniques sont extraits de chaque séquence, puis l'indice de Jaccard est calculé en divisant la taille de l'intersection des k-mers par la taille de leur union.
- Intersection : k-mers communs entre les deux séquences.
- Union : Ensemble total des k-mers présents dans au moins une des deux séquences.

### `if __name__ == "__main__"`

Lorsqu'on compare la même espèce avec elle-même, la similarité de Jaccard sera de par 1 par définition. Sinon, on utilise la fonction jaccard(fileA,fileB,k) pour calcul la similarité. On écrit après également le résultat dans un fichier csv.

## Résultats

### Matrice de Similarité de Jaccard

La matrice de similarité obtenue pour les cinq bactéries est la suivante :

|                      | GCA_000008865.2 | GCA_000005845 | GCA_000069965.1 | GCA_000013265.1 | GCA_030271835.1 |
|----------------------|-----------------|---------------|-----------------|-----------------|-----------------|
| **GCA_000008865.2**   | 1.0             | 0.44601759    | 0.00078665      | 0.30479711      | 0.00077304      |
| **GCA_000005845**     | 0.44601759      | 1.0           | 0.00085698      | 0.33604993      | 0.00085443      |
| **GCA_000069965.1**   | 0.00078665      | 0.00085698    | 1.0             | 0.0008489       | 0.02623378      |
| **GCA_000013265.1**   | 0.30479711      | 0.33604993    | 0.0008489       | 1.0             | 0.00082473      |
| **GCA_030271835.1**   | 0.00077304      | 0.00085443    | 0.02623378      | 0.00082473      | 1.0             |

### Interprétation de la matrice

- **Similitudes élevées** : La similarité la plus forte est observée entre les espèces **GCA_000008865.2** et **GCA_000005845** (Jaccard = 0.446), ainsi que **GCA_000005845** et **GCA_000013265.1** (Jaccard = 0.336). On pourrait faire l'hypothèse qu'elles sont plus proches phylogénétiquement entre elles.
- **Similitudes faibles** : Les valeurs proches de zéro indiquent que les séquences partagent très peu de k-mers. Par exemple, la similarité entre **GCA_000069965.1** et **GCA_000008865.2** est extrêmement faible (Jaccard = 0.00078665), suggérant peu de points communs.

La diagonale contient des 1.0, ce qui est attendu puisque la similarité d'une séquence avec elle-même est toujours maximale.




