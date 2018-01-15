# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:03:52 2018

@author: 3301269
"""

from kmer import *
import matplotlib.pyplot as plt
import sys
import numpy
import os


chaos_directory = "Chaos_archaea_6_100000"
CSV = False
print chaos_directory
liste_genome_directories = sorted(os.listdir("./"+chaos_directory+"/"))

#Recuperer path fichiers fenetres
liste_fichiers = []
Len =[]
liste_dirs = []
for espece in liste_genome_directories :
	print espece
	dossier = os.listdir(chaos_directory+"/"+espece)
	print sorted(dossier)
	print "\n"
	if len(dossier) > 1:
		liste_fichiers += sorted(dossier)[:-1]
		Len.append(len(dossier)-1)
		liste_dirs += [espece for i in xrange(len(dossier)-1)]

def fonction(Len, liste_dirs):
	i = 0
	Index = []
	dirs = []
	j = 0
	while i<len(liste_dirs):
		if Len[j]<4:
			Index.append(range(i,i+Len[j]))
			dirs.append(liste_dirs[i:i+Len[j]])
			print Len[j]
		else:
			indx = sorted(np.random.choice(range(i,i+Len[j]), 4, replace=False))
			Index.append(indx)
			print len(indx)
			for ir in indx:
				dirs.append(liste_dirs[ir])
		i+=Len[j]
		j+=1
	print Index, len(Index)
	return np.concatenate(Index)


def fonction2(Len, liste_dirs):
	i = 0
	Index = []
	dirs = []
	j = 0
	while i<len(liste_dirs):
		#[2,8,41,6,31,38,45,50,24,18]
		# [2,9,52,7,33,40,43,47,26,20]
		# Bact [2,9,49,7,34,42,43,47,26,21]
		if j in [2,8,41,6,31,38,45,50,24,18]:
			if Len[j]<15:
				Index.append(range(i,i+Len[j]))
				dirs.append(liste_dirs[i:i+Len[j]])
				print Len[j]
			else:
				indx = sorted(np.random.choice(range(i,i+Len[j]), 15, replace=False))
				Index.append(indx)
				print len(indx)
				for ir in indx:
					dirs.append(liste_dirs[ir])
		i+=Len[j]
		j+=1
	return Index, dirs
	
def dendrogramme_1_espece_no_clust(chaos_list):
	liste_dendo = []
	for i in range(len(chaos_list)):
		temp = np.ndarray.flatten(chaos_list[i])
		liste_dendo.append(temp)
	return liste_dendo, len(liste_dendo)
		
def dendogramme_sans_cl(liste_chaos_especes, liste_dist_especes):
	"""
	fais un dendrogramme pour les differentes especes
	
	liste_chaos_especes : liste de listes des chaos des fenetres d'une espece
	liste_dist_especes : liste des matrices de distance 2 a 2 des fenetres d'un genome
	"""
	dendo = []
	nom = {}
	i = 0
	cols = ["#000000","#FF0000","#FFFF00","#808000","#00FF00","#008000","#00FFFF","#008080","#0000FF","#000080","#FF00FF","#800080","#C0C0C0","#808080","#800000"]
	r = 0
	for esp in xrange(len(liste_chaos_especes)) : #parcourt les especes
		t, temp = dendrogramme_1_espece_no_clust(liste_chaos_especes[esp])
		dendo += t #Fait une liste de chaos vectorises pour faire le dendogramme
		nom[str(i)] = cols[r]
		for j in range(1,temp):
			nom[str(i+j)] = cols[r]
		i += j+1
		if r + 1 > len(cols):
			print "Pas assez de couleurs! :'("
		r+=1
	print nom
	z = linkage(dendo)
	link_cols = {}
	for i, i12 in enumerate(z[:,:2].astype(int)):
		c1, c2 = (link_cols[x] if x > len(z) else nom["%d"%x]for x in i12)
		link_cols[i+1+len(z)] = c1 
	dendrogram(z,    leaf_rotation=90., leaf_font_size=8.,
				link_color_func=lambda x: link_cols[x])
	print nom
	
	
	
if __name__ == "__main__":
	matrix2 = np.load("matrice_dist_jsdiv_GLOBAL_Chaos_archaea_6_100000.npy")
	matrix = np.load("matrice_dist_jsdiv_FENETRES_Chaos_archaea_6_100000.npy")
	Index = fonction(Len, liste_dirs) #(triche : Len et liste_dirs viennent de generer_distance_matrix_BD.py)
	plt.figure()
	plot_matrice_distance(matrix[Index].T[Index])
	plt.title("Divergence de Jensen-Shannon de\nFenetres genomes vs Fenetres de genomes")
	plt.xlabel("Fenetres de divers genomes")
	plt.ylabel("Fenetres de divers genomes")
	plt.xlim(0, matrix[Index].T[Index].shape[0])	
	plt.ylim(0, matrix[Index].T[Index].shape[0])	
	plt.show()
	plt.figure()
	plot_matrice_distance(matrix2)
	plt.title("Divergence de Jensen-Shannon de\nGenomes vs Genomes")
	plt.xlabel("Signatures de genomes")
	plt.ylabel("Signatures genomes")
	plt.xlim(0, matrix2.shape[0])	
	plt.ylim(0, matrix2.shape[0])	
	plt.show()
#	sous_groupe_fenetres = np.array([  37,   38,   39,   40,   41,   42,   43,   44,   45,   46,   47,
#         48,   49,   50,   51,  100,  101,  102,  103,  104,  105,  106,
#        107,  108,  109,  110,  112,  115,  116,  119,  151,  152,  160,
#        161,  163,  166,  170,  171,  173,  176,  177,  178,  179,  180,
#        181,  344,  345,  346,  349,  351,  353,  354,  355,  356,  360,
#        362,  364,  365,  367,  368,  460,  461,  462,  463,  464,  466,
#        467,  468,  470,  471,  472,  473,  476,  477,  479,  684,  685,
#        686,  687,  688,  689,  690,  691,  692,  693,  694,  695,  696,
#        697,  698,  801,  802,  803,  804,  805,  808,  809,  810,  813,
#        817,  818,  819,  820,  821,  822,  861,  862,  863,  865,  866,
#        867,  868,  869,  870,  871,  872,  873,  874,  875,  876,  929,
#        930,  931,  933,  934,  935,  936,  939,  942,  944,  945,  947,
#        948,  949,  950, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039,
#       1040, 1041, 1042, 1043, 1044, 1045, 1046])
#	sous_groupe_liste_especes = ['Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1','Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome','Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome','Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome','Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA']	
#	sous_matrice = matrix[sous_groupe_fenetres].T[sous_groupe_fenetres]
#	plt.figure()
#	plot_matrice_distance(sous_matrice)
#	plt.show()
	def Recuperer_chaos(liste_fichiers, index, liste_dossiers, dossier_chaos):
		liste_chaos = []
		chaos =[]
		for compt, i in enumerate(index):
			if (compt+1)/15 != (compt+1.)/15.:
				matrix = np.load(chaos_directory + "/" + liste_dossiers[i] +"/"+ liste_fichiers[i])
				chaos.append(matrix)
			else :
				matrix = np.load(chaos_directory + "/" + liste_dossiers[i] +"/"+ liste_fichiers[i])
				chaos.append(matrix)
				liste_chaos.append(chaos)
				chaos = []
		return liste_chaos
	ind, dirs = fonction2(Len, liste_dirs)
	Liste_chaos_esp = Recuperer_chaos(liste_fichiers, np.concatenate(ind), liste_dirs, chaos_directory)
	for i in Liste_chaos_esp:
		print (len(i))
	def Recuperer_distance2a2 (indices, gde_matrice):
		mat = []
		for ind in indices:
			mat.append(gde_matrice[np.array(ind)].T[np.array(ind)])
		return mat
	liste_matrices = Recuperer_distance2a2 (ind, matrix)
#	plt.figure()
#	plot_matrice_distance(liste_matrices[4])
#	plt.show()
	dendrogramme_tous(Liste_chaos_esp,liste_matrices)
	plt.figure()
	dendogramme_sans_cl(Liste_chaos_esp,liste_matrices)
	