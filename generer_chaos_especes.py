# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:40:05 2018

@author: 3301269
"""

# runfile('/users/nfs/Etu9/3301269/Documents/GENOM/generer_chaos_especes.py', args = 'Chaos Genomes', wdir='/users/nfs/Etu9/3301269/Documents/GENOM')

from kmer import *
import sys
import os

if len(sys.argv) != 3 :
	print "Missing output directory argument (1) or genome directory argument (2)"
	sys.exit(1)


directory_chaos = sys.argv[1]
directory_genomes = sys.argv[2]
liste_genomes = os.listdir(directory_genomes) #liste des files dans ce directory

k = 6
fenetres = True
taille = 50000
pas = 50000

for genome in liste_genomes :
	dico = genome_str(directory_genomes + "/" + genome)
	name = "_".join(dico.keys()[0].split(",")[0].split(" ")[1:])
	print name
	os.system("mkdir "+ directory_chaos + "/" + name )
	if fenetres :
		liste_chaos, chaos_genome, kmer_tot_genome, k_mer_list = parcours_genome_fenetre(dico.values()[0], taille, pas, k)
		np.save(directory_chaos + "/" + name + "/chaos_global.npy", chaos_genome)
		for i, chaos in enumerate(liste_chaos):
			np.save(directory_chaos + "/" + name + "/chaos_fenetre_"+ str(i) +".npy", chaos)
	else:
		dico_seq, kmer_tot = kmer_sequence(dico.values()[0], k)
		dico_freq = probabilities_dic (kmer_tot, dico_seq)
		chaos = np.array(chaos_game_representation(dico_freq, k))
		np.save(directory_chaos + "/" + name + "/chaos_global.npy", chaos)

print "Done"