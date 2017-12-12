#Python 2.7
#!/usr/bin/env python 
#-*- coding: utf-8 -*-.

from cgr import chaos_game_representation
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import numpy as np
import pylab
from scipy.stats import chisquare

def genome_str(fichier):
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
	kmer_tot = 0
	for i in xrange(0, len(sequence)-k):
		kmer = sequence[i:i+k]
		if kmer not in dico_kmer.keys():
			dico_kmer[kmer] = 1
			kmer_tot += 1
		else :
			dico_kmer[kmer] += 1
			kmer_tot += 1
	return dico_kmer, kmer_tot



def probabilities_dic (kmer_tot, dico_kmer):
	"""
	Prend un dictionnaire d'occurrences et le tranforme en dictionnaire de frequences

	--------
	input:
	kmer_tot : nombre total de kmers dans la sequence
	dico_kmer : dictionnaire , keys = kmer, value = nombre d'occurences dans la sequence
	
	--------
	output:
	kmer_freq : dictionnaire , keys = kmer, value = frequence de chaque kmer dans la sequence
	"""
	kmer_freq = deepcopy(dico_kmer)
	for key in dico_kmer.keys():
		kmer_freq[key] = kmer_freq[key]/float(kmer_tot)
	return kmer_freq


def plot_CGR (chaos, k):
	"""
	Fonction tiree du main de crg.py qui plot le CRG.

	-------
	input :
	chaos : output de chaos_game_representation() de crg.py
	k : taille k-mer

	-------
	output :
	plot
	"""
	array_size = int(math.sqrt(4**k))
	# Draw the figure with matplotlib

	plt.title('Chaos game representation for ' + str(k) + '-mers' + "\n" ) # gh

	# Assign colors to the values (change cmap for different colors)
	im = plt.imshow(chaos,
		 # To set the plot start at 0
		 extent=[0, array_size, 0, array_size],
		 interpolation='nearest', cmap=cm.hot_r)
	# From http://stackoverflow.com/questions/11776663
	plt.colorbar(im, orientation='vertical')

	# Add black (k) lines to separate the main four regions
	plt.axhline(y=array_size / 2, color='k')
	plt.axvline(x=array_size / 2, color='k')

	# Remove the ticks (unneeded)
	plt.tick_params(
	axis='both', which='both',
	top='off', bottom='off', left='off', right='off',
	labelbottom='off', labelleft='off')

	# Change the ticks to A, T, G and C
	nucleotides = {'A': [0, 1], 'T': [1, 1],
		   'C': [0, 0], 'G': [1, 0]}
	for nuc, i in nucleotides.iteritems():
		plt.annotate(nuc, xy=(i[0], i[1]),
			     xycoords='axes fraction',
			     fontsize=16,
			     xytext=(i[0] * 30 - 20, i[1] * 30 - 23),
			     textcoords='offset points')
   	plt.show()

def euclidean_2chaos (chaos1, chaos2):
	"""
	Calcule la distance euclidienne entre 2 arrays "chaos" qui nous donne la frequence des k-mers.

	-------
	input :
	chaos1 : un array nous donnant la frequence de k-mers sur une sequence
	chaos2 : un array nous donnant la frequence de k-mers sur une autre sequence

	-------
	output :
	float, distance
	"""
	chaos1 = np.ndarray.flatten(chaos1)
	chaos2 = np.ndarray.flatten(chaos2)
	return math.sqrt(np.sum((chaos1-chaos2)**2))


def parcours_genome_fenetre (genome, pas, wind, k) :
	"""
	Parcourt un genome avec un pas de pas et evalue la frequence de k-mers pour une fenetre
	de taille wind.

	-------
	input :
	genome : sequence du genome
	pas : int en pb, pas de parcourt du genome
	wind : int, taille de la fenetre (window)
	k : taille du k-mer

	-------
	output :
	Chaos_list : Liste des arrays chaos (freq k-mer) pour chaque fenetre de genome
	chaos_genome : chaos du genome global (freq k-mer) : signature du genome
	kmer_tot_genome : nombre de k-mers total dans le genome
	k_mer_list : liste du nombre de k-mers dans chaque fenetre
	"""
	Chaos_list = []
	k_mer_list = []
	d_occurr_genome, kmer_tot_genome = kmer_sequence(genome, k)
	d_freq_genome = probabilities_dic (kmer_tot_genome, d_occurr_genome)
	chaos_genome = np.array(chaos_game_representation(d_freq_genome, k))
	for seq in xrange(0, len(genome)-wind, pas) :
		d_occurr, kmer_tot = kmer_sequence(genome[seq:seq+wind], k)
		d_freq = probabilities_dic (kmer_tot, d_occurr)
		chaos = np.array(chaos_game_representation(d_freq, k))
		#plot_CGR (chaos, k)
		Chaos_list.append(chaos)
		k_mer_list.append(kmer_tot)
	return Chaos_list, chaos_genome, kmer_tot_genome, k_mer_list


def matrice_distance2a2 (liste_chaos) :
	"""
	Matrice des distances 2 a 2 euclidienne : fenetres vs fenetres
	
	-------
	input :
	liste_chaos : Liste des arrays chaos (freq k-mer) pour chaque fenetre de genome

	-------
	output :
	dist2a2 : matrice des distances numpy
	"""
	dist2a2 = np.zeros((len(liste_chaos), len(liste_chaos)))
	for i, chaos1 in enumerate(liste_chaos) :
		for j, chaos2 in enumerate(liste_chaos) :
			dist2a2[i,j] = euclidean_2chaos (chaos1, chaos2)
	return dist2a2


def plot_matrice_distance(dist2a2) :	
	"""
	Plot de la matrice des distances 2 a 2
	"""
	pylab.pcolor(dist2a2)
	plt.colorbar()
	plt.title("Carte des distances entre les fenetres du genome vs les fenetres du genome")
	plt.xlabel("Fenetres du genome")
	plt.ylabel("Fenetres du genome")
	plt.show()


def chi_square(chaos_liste, chaos_genome, k_mer_list):
	"""
	Chi2 de conformite des fenetres avec le genome global.

	-------
	input :
	liste_chaos : Liste des arrays chaos (freq k-mer) pour chaque fenetre de genome
	chaos_genome : chaos du genome global (freq k-mer) : signature du genome
	k_mer_list : liste du nombre de k-mers dans chaque fenetre
	"""
	pval_liste = []
	for i, chaos in enumerate(chaos_liste): # /!\ Effectifs superieurs a 5 pour le theorique et l'observe
		if np.any(chaos * k_mer_list[i] < 5) :
			print "Chi-square ERROR : Effectifs observes inferieurs a 5!"
			break
		elif np.any(chaos_genome * k_mer_list[i] < 5) :
			print "Chi-square ERROR : Effectifs theoriques inferieurs a 5!"
			break
		pval_liste.append(chisquare(chaos * k_mer_list[i], (k_mer_list[i] * chaos_genome).astype(int), axis = None))
	return pval_liste

def distance_a_ref (liste_chaos, chaos_genome):
	"""
	Distances 2 a 2 euclidienne : fenetres vs genome global
	
	-------
	input :
	liste_chaos : Liste des arrays chaos (freq k-mer) pour chaque fenetre de genome
	chaos_genome : chaos du genome global (freq k-mer) : signature du genome

	-------
	output :
	vecteur des distances numpy
	"""
	dist = []
	for i, chaos1 in enumerate(liste_chaos) :
		dist.append(euclidean_2chaos (chaos1, chaos_genome))
	return np.array(dist)

dico = genome_str("GCA_000017165.1_ASM1716v1_genomic.fna")
dico2 = genome_str("GCA_000762265.1_ASM76226v1_genomic.fna")
dico3 = genome_str("GCA_000953115.1_DSM1535_genomic.fna")

print len(dico.values()[0]), len(dico2.values()[0]), len(dico3.values()[0])
#2.494.510 2.449.987 2.478.074

a = "CCTCACCAGCGGAAAGTTTAAATATGGATACCATACAATTTTTAGTACCTATTGCAATCTGCGGTGGATCCGCTCACATTGTATGCCCTGATACGATGTGGTCTGTGGTTTTACTTGCCCAGTCTGCATTGTGGAAATATTTATAAATAGATCCGGACAGATATTAATAGATGAATAGAGTAGATTTGTCCATATTTATCCCGGATTCACTGACGGCTGAGACAGGGGATCTCAAAATAAAGACCTACAAGGTGGGTCTCATTGCACGGGCCGCTTCGATATTCGGGGTTAAGCGTATAGTGATCTATCACGATGATGCAGATGGAGAGGCAAGGTTTATTAGGGATATCCTGACTTATATGGATACACCTCAATACCTTCGCAGGAAGGTTTTCCCGATAATGAGGGAGTTGAAACATGTGGGGATACTCCCACCTCTGAGAACTCCCCATCACCCAACCGGAAAACCCGTTACTGGTGAATACAGACAGGGACTGACAGTTAAAAGGGTAAAGAAAGGAACTCTTGTGGATATTGGCGCAGATAAACTTGCACTGTGCAGGGAAAAACTGACAGTGAATAGGATAATGAGTTTCAGGGTTGTCAGGTTGGGTAAGGAAATACTGATAGAGCCCGATGAACCAGACGATAGATACTGGGGATACGAGGTACTGGATACCCGGAGGAACCTCGCAGAGAGCCTTAAAACATTAGGTGCCGATGTTGTCGTGGCAACATCCAGGAAAGCTTCGCCCATTACTTCTATTCTGGATGAAGTAAAAACGAGGATGAGGGGGGCC"

dico_seq, kmer_tot = kmer_sequence(dico2.values()[0],2)

dico_freq = probabilities_dic (kmer_tot, dico_seq)

chaos2 = np.array(chaos_game_representation(dico_freq, 2))
plot_CGR (chaos2, 2)

liste_chaos, chaos_genome, kmer_tot_genome, k_mer_list = parcours_genome_fenetre(dico.values()[0], 100000, 100000, 2)
print chi_square(liste_chaos, chaos_genome, k_mer_list)
print chi_square(liste_chaos, chaos2, k_mer_list)
print len(liste_chaos)
dist2a2 = matrice_distance2a2(liste_chaos)
print dist2a2, dist2a2.shape
plot_matrice_distance(dist2a2)

distancesAref = distance_a_ref (liste_chaos, chaos_genome)
plt.bar(range(0, len(distancesAref)), distancesAref)
plt.ylim(0,0.1)
plt.show()

print distancesAref


