import numpy as np
import heapq
import random
from TP.kmers import stream_kmers, kmer2str
import hashlib
def filter_smallest(seq,k,s):
	lst=[np.inf for _ in range(s)]
	max=np.inf
	for kmer in stream_kmers(seq,k):
		if kmer<max_element:
			idx=lst.index(max_element)
			lst[idx]=kmer
			max_element=max(lst)
	return lst


#def filter_smallest(seq, k, s):
	"""Garde les s plus petits k-mers d'une séquence."""
#	# On initialise une liste vide pour le heap
#	smallest_kmers = []
	# On itère sur les k-mers de la séquence
#	for kmer in stream_kmers(seq, k):
		# On inverse kmer pour utiliser heapq comme un max-heap 
		# On stocke les éléments sous forme négative pour simuler un max-heap
#		if len(smallest_kmers) < s:
#			heapq.heappush(smallest_kmers, -kmer)  # Ajoute le k-mer au tas
#		else:
			# Si le tas a déjà s éléments, on remplace le plus grand si le nouveau k-mer est plus petit
#			if -kmer > smallest_kmers[0]:
#				heapq.heapreplace(smallest_kmers, -kmer)

	# On renverse les éléments à la fin pour les remettre dans leur valeur originale
#	return [-k for k in smallest_kmers]
    
def hyperMinSketch(seq1, seq2):

	i, j = 0, 0
	intersection = 0
	union = 0
    
	while i < len(seq1) and j < len(seq2):
		if seq1[i] == seq2[j]:
			intersection += 1
			union += 1
			i += 1
			j += 1
		elif seq1[i] < seq2[j]:
			union += 1
			i += 1
		else:
			union += 1
			j += 1
    
	while i < len(seq1):
		union += 1
		i += 1
    
	while j < len(seq2):
		union += 1
		j += 1
    
	return intersection, union
    
    
def xorshift(val):
	"""Générateur de nombres pseudo-aléatoires basé sur xorshift"""
	val ^= val << 13
	val &= 0xFFFFFFFFFFFFFFFF
	val ^= val >> 17
	val &= 0xFFFFFFFFFFFFFFFF
	val ^= val << 5
	return val




def hash_kmer(kmer, i):
	"""Retourne une valeur de hachage pour le k-mer combiné à l'indice i."""
	hash_object = hashlib.sha1(f"{kmer}_{i}".encode())
	return int(hash_object.hexdigest(), 16)
	
	
def partition_min_hash(kmer_list, num_partitions, num_hashes):
	partition_size = 1
	partitions = [kmer_list[i:i + partition_size] for i in range(0, len(kmer_list), partition_size)]
	partitions = partitions[:num_partitions]
	min_hash_signatures = []
	for partition in partitions:
		min_hashes = [float('inf')] * num_hashes
		for kmer in partition:
			for i in range(num_hashes):
				#hash_value = xorshift(kmer ^ i)
				hash_value=hash_kmer(kmer,i)
				#print(f"k-mer: {kmer}, Hash {i}: {hash_value}")
				if hash_value < min_hashes[i]:
					min_hashes[i] = hash_value
		min_hash_signatures.append(min_hashes)
    
	return min_hash_signatures[:num_partitions]


    
