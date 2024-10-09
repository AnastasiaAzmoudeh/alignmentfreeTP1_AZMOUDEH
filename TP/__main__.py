from TP.loading import load_directory,load_fasta
from TP.kmers import stream_kmers, kmer2str,enumerate_kmers
import numpy as np
import csv


def jaccard(fileA, fileB, k):
    j = 0
    # --- To complete ---
    
    """Calcule la similarité de Jaccard entre les k-mers de deux fichiers FASTA."""
    
    # Charger les séquences des fichiers
    
    sequencesA = load_fasta(fileA)
    sequencesB = load_fasta(fileB)
    
    # Concaténer toutes les séquences de l'échantillon en une seule chaîne
    seqA = ''.join(sequencesA)  
    seqB = ''.join(sequencesB)
    
    # Énumérer les k-mers des deux séquences
    kmersA = set(enumerate_kmers(seqA, k))  # On utilise un set pour avoir les k-mers uniques
    kmersB = set(enumerate_kmers(seqB, k))  # Idem
    
    # Initialisation des compteurs
    intersection_count = 0
    union_count = 0
    
    # Calcul de l'intersection (les k-mers communs)
    for kmer in kmersA:
        if kmer in kmersB:
            intersection_count += 1  # On incrémente le compteur pour l'intersection
    
    # Calcul de l'union (les k-mers uniques entre les deux séquences)
    union_set = kmersA.union(kmersB)  # Unir les deux ensembles de k-mers
    union_count = len(union_set)  # On compte tous les k-mers dans l'union
    
    # Calcul de la similarité de Jaccard
    if union_count == 0:
        return 0.0  # Éviter la division par zéro si l'union est vide
    j=intersection_count / union_count  # Similarité de Jaccard
    return j



if __name__ == "__main__":
   
	print("Computation of Jaccard similarity between files")

	# Load all the files in a dictionary
	files = load_directory("data")
	k = 21
    
	print("Computing Jaccard similarity for all pairs of samples")
	filenames = list(files.keys())
	num_files = len(filenames)
	
	print("Loaded files:")
	for filename, sequences in files.items():
		print(f"Filename: {filename}, Number of sequences: {len(sequences)}")
		print(sequences)

	# Créer une matrice de similarité de Jaccard (initialisée avec des zéros)
	jaccard_matrix = np.zeros((num_files, num_files))
	for i in range(len(files)):
		for j in range(i, len(files)):
            
			# --- Complete here ---
            
			if filenames[i] == filenames[j]:
				# La similarité d'un fichier avec lui-même est toujours 1
				jaccard_matrix[i, j] = 1.0
				jaccard_matrix[j, i] = 1.0
			else:
                
				# On calcule l'indice de Jaccard
				jac= jaccard(files[filenames[i]], files[filenames[j]], k)
				print("jac =",jac)
				print(filenames[i],filenames[j],jac)

				# On Remplit la matrice de Jaccard (symétrique)
				jaccard_matrix[i, j] = jac
				jaccard_matrix[j, i] = jac  # Symétrique

	# Affichage de la matrice de Jaccard
	print("Jaccard similarity matrix:")
	print(jaccard_matrix)

    # Sauvegarder la matrice dans un fichier CSV
	output_file = "jaccard_matrix.csv"
	with open(output_file, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile)

		# Écrire l'entête (noms des échantillons)
		writer.writerow([''] + filenames)

		# Écrire les lignes de la matrice
		for i in range(num_files):
			writer.writerow([filenames[i]] + list(jaccard_matrix[i, :]))

		print(f"Jaccard similarity matrix saved to {output_file}")

            
