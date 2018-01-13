# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 11:34:57 2018

@author: 3301269
"""

#runfile('/Vrac/GENOM_adele/generer_distance_matrix_par_espece.py', args='Chaos_bacteria_6_100000', wdir='/Vrac/GENOM_adele')

from kmer import *
import sys
import os

if len(sys.argv) != 2 :
	print 'Missing genome directory argument (without "/") (1)'
	sys.exit(1)


chaos_directory = sys.argv[1]
liste_genome_directories = sorted(os.listdir(chaos_directory)) #liste des dossiers d'especes dans ce directory par ordre alphabetique
fenetres = True

#Dossier chaos de la forme Chaos/nom_espece/chaos_fenetres et chaos_global

#Recuperer matrices	
def recuperer_matrices(chaos_directory, liste_fichiers, liste_dirs):
	liste_chaos = []
	for i, fichier in enumerate(liste_fichiers) :
		matrix = np.load(chaos_directory + "/" + liste_dirs[i] +"/"+ fichier)
		liste_chaos.append(matrix)
	return liste_chaos

os.system("mkdir Matrice_"+ chaos_directory)

#Recuperer path fichiers fenetres
for espece in liste_genome_directories :
	liste_fichiers = []
	liste_dirs = []
	dossier = os.listdir(chaos_directory+"/"+espece)
	if len(dossier) > 1:
		liste_fichiers += sorted(dossier)[:-1]
		liste_dirs += [espece for i in xrange(len(dossier)-1)]
	#Matrice distance fenetres
	liste_chaos_fenetres = recuperer_matrices(chaos_directory, liste_fichiers, liste_dirs)
	matrice_eucli_fen = matrice_distance2a2 (liste_chaos_fenetres)
	np.save("Matrice_"+ chaos_directory + "/matrice_dist_eucli_FENETRES_"+chaos_directory+"_" + espece + ".npy", matrice_eucli_fen)
	matrice_js_fen = matrice_distance2a2_jsdiv (liste_chaos_fenetres)
	np.save("Matrice_"+ chaos_directory + "/matrice_dist_jsdiv_FENETRES_"+chaos_directory+"_" + espece + ".npy", matrice_js_fen)