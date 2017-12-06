#Python 2.7
#!/usr/bin/env python 
#-*- coding: utf-8 -*-.

def sequence_str(fichier):
	"""
	transforme une fichier fasta d'un genome en un str

	---------
	input :
	fichier : fichier fasta (.fna)

	---------
	output:
	dico : dictionnaire , key = nom de l espece, value = sequence
	"""
	f = open(fichier,"r")
	ligne = f.readline()
	dico = {}
	while ligne != "":
		if ligne[0] == ">":
			nom = ligne[1:].rstrip("\n")
			dico[nom] = ""
			ligne = f.readline()	
		else :
			dico[nom]+=ligne.rstrip("\n")	
			ligne = f.readline()		
	f.close()
	return dico

#dico = sequence_str("GCA_002356395.1_ASM235639v1_genomic.fna")
a = "CCTCACCAGCGGAAAGTTTAAATATGGATACCATACAATTTTTAGTACCTATTGCAATCTGCGGTGGATCCGCTCACATTGTATGCCCTGATACGATGTGGTCTGTGGTTTTACTTGCCCAGTCTGCATTGTGGAAATATTTATAAATAGATCCGGACAGATATTAATAGATGAATAGAGTAGATTTGTCCATATTTATCCCGGATTCACTGACGGCTGAGACAGGGGATCTCAAAATAAAGACCTACAAGGTGGGTCTCATTGCACGGGCCGCTTCGATATTCGGGGTTAAGCGTATAGTGATCTATCACGATGATGCAGATGGAGAGGCAAGGTTTATTAGGGATATCCTGACTTATATGGATACACCTCAATACCTTCGCAGGAAGGTTTTCCCGATAATGAGGGAGTTGAAACATGTGGGGATACTCCCACCTCTGAGAACTCCCCATCACCCAACCGGAAAACCCGTTACTGGTGAATACAGACAGGGACTGACAGTTAAAAGGGTAAAGAAAGGAACTCTTGTGGATATTGGCGCAGATAAACTTGCACTGTGCAGGGAAAAACTGACAGTGAATAGGATAATGAGTTTCAGGGTTGTCAGGTTGGGTAAGGAAATACTGATAGAGCCCGATGAACCAGACGATAGATACTGGGGATACGAGGTACTGGATACCCGGAGGAACCTCGCAGAGAGCCTTAAAACATTAGGTGCCGATGTTGTCGTGGCAACATCCAGGAAAGCTTCGCCCATTACTTCTATTCTGGATGAAGTAAAAACGAGGATGAGGGGGGCC"

def kmer_sequence(sequence,k):
	"""
	prend un genome et cree un dictionnaire d'occurence de k mer

	---------
	input :
	sequence : str, sequence nt
	k : int, taille du kmer

	---------
	output:
	dico_kmer : dictionnaire , keys = kmer, value = nombre d'occurences dans la sequence
	"""
	dico_kmer = {}
	for i in xrange(0,len(sequence)-k):
		kmer = sequence[i:i+k]
		if kmer not in dico_kmer.keys():
			dico_kmer[kmer] = 1
		else :
			dico_kmer[kmer] += 1
	return dico_kmer

print kmer_sequence(a,2)




