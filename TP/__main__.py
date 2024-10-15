from TP.loading import load_directory,load_fasta
from TP.kmers import stream_kmers, kmer2str
from TP.hash import filter_smallest,hyperMinSketch,xorshift,partition_min_hash


def jaccard_from_sorted_lists(lstA, lstB):
    idxA = 0
    idxB = 0

    intersection = 0
    union = 0

    while idxA < len(lstA) and idxB < len(lstB):
        union += 1
        if lstA[idxA] == lstB[idxB]:
            intersection += 1
            idxA += 1
            idxB += 1
        elif lstA[idxA] < lstB[idxB]:
            idxA += 1
        else:
            idxB += 1

    union += len(lstA) - idxA
    union += len(lstB) - idxB

    return intersection / union

def jaccard_hyperMinSketch_from_files(file1, file2, k,s):

    
    seq1 = load_fasta(file1)
    seq2 = load_fasta(file2)
    
    # On applique la technique de filtrage pour ne garder que les plus petits k-mers
    seq1 = filter_smallest(seq1,k,s)
    seq2 = filter_smallest(seq2,k,s)
    
    for s1 in seq1:
        kmer_list1.extend(filter_smallest(s1, k, s))  # Par exemple, garder les 100 plus petits k-mers
    for s2 in seq2:
        kmer_list2.extend(filter_smallest(s2, k, s))
    # On utilise hyperMinSketch pour calculer l'intersection et l'union des k-mers
    intersection, union = hyperMinSketch(seq1, seq2)
    
    # Calcul de la similarité de Jaccard
    if union == 0:
        return 0.0  # Si l'union est vide, on retourne 0 pour éviter la division par zéro.
    
    jaccard_similarity = intersection / union
    return jaccard_similarity
    
def jaccard_partition_minhash(seq1, seq2, num_partitions, num_hashes,s,k):
	#seq1 = filter_smallest(seq1,k,s)
	#seq2 = filter_smallest(seq2,k,s)
	minhash_seq1 = partition_min_hash(seq1, num_partitions, num_hashes)
	minhash_seq2 = partition_min_hash(seq2, num_partitions, num_hashes)
	intersection = 0
	union = 0
    
	# on compare les partitions une par une
	for i in range(num_partitions):
		for h in range(num_hashes):
			union += 1
			if minhash_seq1[i][h] == minhash_seq2[i][h]:
				intersection += 1
    
    
	return intersection / union

if __name__ == "__main__":
	print("Computation of Jaccard similarity between files")

	# Load all the files in a dictionary
	print("Loading files")
	files = load_directory("data")
	k = 21
	s=10000
	num_partitions=10000
	num_hashes=10000
	
	
    
	filenames = list(files.keys())
    
	# Create all the kmer lists (can be expensive in memory)
	print("Computing all kmer vectors")
	kmer_lists = {}
	for filename in filenames:
		kmer_lists[filename] = []
		# Enumerate all the sequences from a fasta
		for seq in files[filename]:
			kmer_lists[filename].extend(stream_kmers(seq, k))
		# Sort the kmer lists to speed up the comparison
		kmer_lists[filename].sort()

	#print("Computing Jaccard similarity for all pairs of samples")
	#for i in range(len(files)):
		#for j in range(i+1, len(files)):
			#jaccard = jaccard_from_sorted_lists(kmer_lists[filenames[i]], kmer_lists[filenames[j]])
            
	#        print(filenames[i], filenames[j], jaccard)
	
	#print("Testing Jaccard similarity with HyperMinSketch")
	#for i in range(len(files)):
		#for j in range(i+1, len(files)):
			#file1 = filenames[i]
			#file2 = filenames[j]
			#jaccard_hyperminsketch = jaccard_hyperMinSketch_from_files(file1, file2,k,1000)
			#print(f"{filenames[i]} vs {filenames[j]} - HyperMinSketch Jaccard: {jaccard_hyperminsketch}")
            
                   
	print("Testing Jaccard similarity with Partition MinHash")
	for i in range(len(files)):
		for j in range(i+1, len(files)):
			seq1 = kmer_lists[filenames[i]]
			seq2 = kmer_lists[filenames[j]]
			jaccard_minhash = jaccard_partition_minhash(seq1, seq2, num_partitions, num_hashes,s,k)
			print(f"{filenames[i]} vs {filenames[j]} - Partition MinHash Jaccard: {jaccard_minhash}")
