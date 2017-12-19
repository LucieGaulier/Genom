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
from sklearn.decomposition import PCA

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


def tri_signature(dist2a2, seuil):
	"""
	Trouve les fenetres contenant au moins la moitie de ces distances superieur a un seuil

	-----
	input :
	seuil -float, pourcentage (0.75 pour 75%) de l'etendue des distances (max et min) au 		dessus de laquelle on considere que la fenetre est eloignee
	dist2a2 - matrice des distances 2 a 2 numpy entre les fenetres du genome

	-----
	ouput :
	liste des index des fenetres
	"""
	maxi = np.max(dist2a2)
	dist2a2[dist2a2 ==0] = np.inf
	mini = np.min(dist2a2)
	dist2a2[dist2a2 == np.inf] = 0
	liste = []
	chiffre_seuil = (maxi-mini)*seuil+mini
	for i in range(0,len(dist2a2)):
		n = 0.
		for j in range(0,len(dist2a2[i])):
			#print dist2a2[i][j],'ok'
			if dist2a2[i][j] >= chiffre_seuil :
				n+=1
		temp = n/len(dist2a2[i])
		if temp >= 0.5 :
			liste.append(i)
	return liste,mini

def trouver_transferts_moyenne(dist2a2, seuil):
	"""
	Trouve les fenetres contenant des moyennes de distance aux autres fenetres superieures
	a un seuil.

	-------
	input :
	dist2a2 - matrice des distances 2 a 2 numpy entre les fenetres du genome
	seuil - float, pourcentage (0.75 pour 75%) de l'etendue des distances (max et min) au dessus de laquelle on considere que la fenetre est eloignee

	------
	output :
	liste des index des fenetres
	"""
	dist2a2[dist2a2 ==0] = np.inf
	mini = np.min(dist2a2)
	dist2a2[dist2a2 == np.inf] = 0
	seuil = (np.max(dist2a2)-mini) * seuil + mini
	moyennes = np.mean(dist2a2, axis = 0)
	return np.where(moyennes >= seuil)

def clustering_transferts_horizontaux(liste_index, dist2a2,seuil,mini):
	"""
	fait du clustering entre les fenetres potentiellement issues d'un transfert : met ensemble les fenetres inferieures a un seuil, donc hypothetiquement issus du meme transfert

	-------
	inputs :
	liste_index : liste des index des fenetres "de transfert"
	dist2a2 : matrice des distances 2 a 2 numpy entre toutes les fenetres du genome
	seuil : float, pourcentage (0.75 pour 75%) de l'etendue des distances (max et min) au dessus de laquelle on considere que la fenetre est eloignee
	mini : minimum de la matrice (a part 0)

	-------
	output :
	clusters : liste de set, contenant les indices des fenetres clusterisees
	"""
	maxi = np.max(dist2a2)
	liste_index = np.array(liste_index)
	chiffre_seuil = (maxi-mini)*seuil+mini
	liste = dist2a2[liste_index].T[liste_index].T
	clusters = []
	for i, fen in enumerate(liste): #Parcourt les distances de la fenetre vs toutes les autres 
		proches_i = liste_index[fen < chiffre_seuil] #indices des fenetres proches (inf seuil) de celle que l'on parcourt
		in_cluster = False
		for cl in xrange(len(clusters)):
			if np.any(np.isin(proches_i, list(clusters[cl]))) :
				clusters[cl] = clusters[cl].union(proches_i)
				in_cluster = True
		if not in_cluster:
			clusters.append(set(proches_i))
	return clusters
			
def pca(X,Y,n):
    pca = PCA(n_components=n)
    pca.fit(X)
    X_pca = pca.transform(X)
    X_back = pca.inverse_transform(X_pca)
#    plt.figure()
#    plt.title("projection by PCA")
#    plt.scatter(X_pca[:, 0], X_pca[:, 1],marker='o', c=Y,s=25, edgecolor='k')
#    plt.show()
    return X_pca
    
"""
#chaos = freq
#on veut trouver les chaos qui sont les plus differents de tout les autres
#on trouve les transferts : different de plus de la moitie des autres fenetre
#si on en trouve plusieurs > on compare les transferts entre eux pour savoir si ils viennent du meme transfert ou non

Recuperer signature des transferts et du genome.
genome > le plus grand ? le majoraitaire, on part du principe que si transfert c'est "petite" patate

cherche signature transfert dans les autres génome (BD)

si plsr fois la meme espece > tri ?

parfois N dans sequence

Revoir un peu plus pour lhistoire de signature globale = signature parties
"""

	
#dico = genome_str("GCA_000017165.1_ASM1716v1_genomic.fna")
dico2 = genome_str("GCA_000762265.1_ASM76226v1_genomic.fna")
dico3 = genome_str("GCA_000953115.1_DSM1535_genomic.fna")
dico_cool = genome_str("GCA_000016525.1_ASM1652v1_genomic.fna")
dico_chelou = genome_str("GCA_001889405.1_ASM188940v1_genomic.fna")

print len(dico_cool.values()[0]), len(dico2.values()[0]), len(dico3.values()[0])
#2.494.510 2.449.987 2.478.074
k=5
a = "CCTCACCAGCGGAAAGTTTAAATATGGATACCATACAATTTTTAGTACCTATTGCAATCTGCGGTGGATCCGCTCACATTGTATGCCCTGATACGATGTGGTCTGTGGTTTTACTTGCCCAGTCTGCATTGTGGAAATATTTATAAATAGATCCGGACAGATATTAATAGATGAATAGAGTAGATTTGTCCATATTTATCCCGGATTCACTGACGGCTGAGACAGGGGATCTCAAAATAAAGACCTACAAGGTGGGTCTCATTGCACGGGCCGCTTCGATATTCGGGGTTAAGCGTATAGTGATCTATCACGATGATGCAGATGGAGAGGCAAGGTTTATTAGGGATATCCTGACTTATATGGATACACCTCAATACCTTCGCAGGAAGGTTTTCCCGATAATGAGGGAGTTGAAACATGTGGGGATACTCCCACCTCTGAGAACTCCCCATCACCCAACCGGAAAACCCGTTACTGGTGAATACAGACAGGGACTGACAGTTAAAAGGGTAAAGAAAGGAACTCTTGTGGATATTGGCGCAGATAAACTTGCACTGTGCAGGGAAAAACTGACAGTGAATAGGATAATGAGTTTCAGGGTTGTCAGGTTGGGTAAGGAAATACTGATAGAGCCCGATGAACCAGACGATAGATACTGGGGATACGAGGTACTGGATACCCGGAGGAACCTCGCAGAGAGCCTTAAAACATTAGGTGCCGATGTTGTCGTGGCAACATCCAGGAAAGCTTCGCCCATTACTTCTATTCTGGATGAAGTAAAAACGAGGATGAGGGGGGCC"

dico_seq, kmer_tot = kmer_sequence(dico2.values()[0],k)

dico_freq = probabilities_dic (kmer_tot, dico_seq)

chaos2 = np.array(chaos_game_representation(dico_freq,k))
plot_CGR (chaos2, k)

liste_chaos, chaos_genome, kmer_tot_genome, k_mer_list = parcours_genome_fenetre(dico_chelou.values()[0], 50000,50000, k)
#print chi_square(liste_chaos, chaos_genome, k_mer_list) #on veut des p values grandes (non rejet de H0
#print chi_square(liste_chaos, chaos2, k_mer_list) #on veut des p values petite (rejet de H0)

#print liste_chaos
dist2a2 = matrice_distance2a2(liste_chaos)
#print dist2a2, dist2a2.shape

plt.figure()
plot_matrice_distance(dist2a2)
plt.show()
"""
distancesAref = distance_a_ref (liste_chaos, chaos_genome)
gen = euclidean_2chaos(chaos_genome, chaos2)

plt.figure()
plt.bar(range(0, len(distancesAref)), distancesAref)
plt.plot((0, len(distancesAref)),(gen,gen),'k-')
plt.ylim(0,0.1)
plt.show()
"""
#print distancesAref
#print gen	

#print trouver_transferts_moyenne(dist2a2, 0.60)

index,m = tri_signature(dist2a2, 0.45)
print index
b = clustering_transferts_horizontaux(index,dist2a2,0.45,m)
print "clusters", b

liste_y = [0 for i in range(len(dist2a2))]
for j in range(1,len(b)+1):
    for k in range(len(b[j-1])):
        liste_y[list(b[j-1])[k]] = j
        
#liste_y[-1] = 3

z = pca(dist2a2,liste_y,2)
plt.figure()
plt.plot(z)
plt.show()
plt.figure()
plt.scatter(z[:,0],z[:,1],marker='o',s=25,c=liste_y)
plt.show()
