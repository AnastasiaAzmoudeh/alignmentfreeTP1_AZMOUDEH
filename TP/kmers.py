
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def encode_nucl(lettre):
	"""Encode un nucléotide sous forme d'entier (A=0, C=1, T=2, G=3)."""
	encoding = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
	return encoding[lettre]
	
def encode_kmer(seq, k):
	"""Encode un k-mer sous forme d'entier."""
	kmer = 0
	for letter in seq[:k]:
		kmer = (kmer << 2) | encode_nucl(letter)  # Décale et ajoute le nucléotide
	return kmer

def enumerate_kmers(seq, k):
	"""Génère tous les k-mers d'une séquence sous forme de chaînes de caractères."""
	mask = (1 << (2 * k)) - 1  # Masque pour conserver les k derniers nucléotides encodés
	kmer = encode_kmer(seq[:k], k)  # Encode le premier k-mer
	yield kmer2str(kmer, k)  # Génère la version chaîne de caractères du premier k-mer
	
	for i in range(1, len(seq) - k + 1):
		kmer = ((kmer << 2) & mask) | encode_nucl(seq[i + k - 1])  # Décale et ajoute le nouveau nucléotide
		yield kmer2str(kmer, k)  # Génère la version chaîne de caractères du k-mer

	
def reverse_complement(kmer):
	"""Génére le complément inverse d'un k-mer"""
	complement={'A':'T','C':'G','G':'C','T':'A'}
	return ''.join(complement[nuc] for nuc in reversed(kmer))
    
    
def stream_kmers(text, k):
	"""Renvoie les k-mer canoniques"""
    # --- To complete ---
	for kmer in enumerate_kmers(text, k):
		rev_comp = reverse_complement(kmer)
		canonical = min(kmer, rev_comp)
		yield canonical
    
    #pass
