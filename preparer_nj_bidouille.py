# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:03:52 2018

@author: 3301269
"""

from kmer import *
import sys
import numpy
import os

CSV = False

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
	matrix2 = np.load("matrice_dist_eucli_GLOBAL_Chaos_archaea_6_100000.npy")
	matrix = np.load("matrice_dist_eucli_FENETRES_Chaos_archaea_6_100000.npy")
	Index = fonction(Len, liste_dirs) #(triche : Len et liste_dirs viennent de generer_distance_matrix_BD.py)
	plt.figure()
	plot_matrice_distance(matrix[Index].T[Index])
	plt.title("Distances euclidiennes de\nFenetres genomes vs Fenetres de genomes")
	plt.xlabel("Fenetres de divers genomes")
	plt.ylabel("Fenetres de divers genomes")
	plt.xlim(0,198)
	plt.ylim(0,198)
	plt.show()
	plt.figure()
	plot_matrice_distance(matrix2)
	plt.title("Distances euclidiennes de\nGenomes vs Genomes")
	plt.xlabel("Signatures de genomes")
	plt.ylabel("Signatures genomes")
	plt.xlim(0,53)
	plt.ylim(0,53)
	plt.show()
	sous_groupe_fenetres = np.array([  37,   38,   39,   40,   41,   42,   43,   44,   45,   46,   47,
         48,   49,   50,   51,  100,  101,  102,  103,  104,  105,  106,
        107,  108,  109,  110,  112,  115,  116,  119,  151,  152,  160,
        161,  163,  166,  170,  171,  173,  176,  177,  178,  179,  180,
        181,  344,  345,  346,  349,  351,  353,  354,  355,  356,  360,
        362,  364,  365,  367,  368,  460,  461,  462,  463,  464,  466,
        467,  468,  470,  471,  472,  473,  476,  477,  479,  684,  685,
        686,  687,  688,  689,  690,  691,  692,  693,  694,  695,  696,
        697,  698,  801,  802,  803,  804,  805,  808,  809,  810,  813,
        817,  818,  819,  820,  821,  822,  861,  862,  863,  865,  866,
        867,  868,  869,  870,  871,  872,  873,  874,  875,  876,  929,
        930,  931,  933,  934,  935,  936,  939,  942,  944,  945,  947,
        948,  949,  950, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039,
       1040, 1041, 1042, 1043, 1044, 1045, 1046])
	sous_groupe_liste_especes = ['Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Archaeoglobus_profundus_DSM_5631', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1','Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Halobacterium_sp._NRC-1', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Haloquadratum_walsbyi_DSM_16790_complete_genome', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanococcoides_burtonii_DSM_6242', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome','Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanohalophilus_halophilus_strain_Z-7982_chromosome', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Methanosphaera_stadtmanae_DSM_3091', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrobaculum_aerophilum_str._IM2', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome','Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome','Pyrococcus_abyssi_complete_genome', 'Pyrococcus_abyssi_complete_genome', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Sulfolobus_acidocaldarius_DSM_639', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA', 'Thermoplasma_volcanium_GSS1_DNA']	
	sous_matrice = matrix[sous_groupe_fenetres].T[sous_groupe_fenetres]
	plot_matrice_distance(sous_matrice)
	plt.show()
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
	Liste_chaos_esp = Recuperer_chaos(liste_fichiers, sous_groupe_fenetres, liste_dirs, "Chaos_archaea_6_100000")
	for i in Liste_chaos_esp:
		print (len(i))
	ind, dirs = fonction2(Len, liste_dirs)
	def Recuperer_distance2a2 (indices, gde_matrice):
		mat = []
		for ind in indices:
			mat.append(gde_matrice[np.array(ind)].T[np.array(ind)])
		return mat
	liste_matrices =  Recuperer_distance2a2 (ind, matrix)
	plt.figure()
	plot_matrice_distance(liste_matrices[4])
	plt.show()
	dendrogramme_tous(Liste_chaos_esp,liste_matrices)
	plt.figure()
	dendogramme_sans_cl(Liste_chaos_esp,liste_matrices)
	if CSV :
		names_bac = ["Acholeplasma_laidlawii_PG-8A","Acidobacterium_capsulatum_ATCC_51196","Akkermansia_muciniphila_ATCC_BAA-835","Alicyclobacillus_acidocaldarius_subsp._acidocaldarius_DSM_446_plasmid_pAACI02","Aquifex_aeolicus_VF5","Bacillus_cereus_Q1_plasmid_pBc53","Bacillus_pseudofirmus_OF4","Bacteroides_fragilis_YCH46_DNA","Bdellovibrio_bacteriovorus_complete_genome","Bordetella_pertussis_strain_Tohama_I","Borrelia_burgdorferi_B31_plasmid_cp32-6","Campylobacter_jejuni_subsp._jejuni_81-176_plasmid_pTet","Candidatus_Amoebophilus_asiaticus_5a2","Candidatus_Cloacamonas_acidaminovorans_str._Evry_provisional_genome_sequence_from_WWE1_candidate_division","Candidatus_Endomicrobium_trichonymphae_plasmid_pRSTT1_DNA","Candidatus_Solibacter_usitatus_Ellin6076_chromosome","Carboxydothermus_hydrogenoformans_Z-2901","Chlamydia_trachomatis_strain_L2_434_Bu_complete_genome","Chlorobium_chlorochromatii_CaD3","Chloroflexus_aurantiacus_J-10-fl","Clostridium_acetobutylicum_ATCC_824_megaplasmid_pSOL1","Corynebacterium_glutamicum_ATCC_13032","Corynebacterium_glutamicum_ATCC_13032_DNA","Coxiella_burnetii_RSA_493_plasmid_pQpH1","Cyanothece_sp._ATCC_51142_plasmid_C","Dehalococcoides_mccartyi_195_chromosome","Deinococcus_radiodurans_R1_chromosome_2","Deinococcus_radiodurans_R1_plasmid_MP1","Dictyoglomus_thermophilum_H-6-12","Elusimicrobium_minutum_Pei191","Fibrobacter_succinogenes_subsp._succinogenes_S85","Flavobacterium_psychrophilum_JIP02_86_complete_genome","Fusobacterium_nucleatum_subsp._nucleatum_ATCC_25586","Gemmatimonas_aurantiaca_T-27_DNA","Gloeobacter_violaceus_PCC_7421_DNA","Leptospira_interrogans_serovar_Lai_str._56601_chromosome_I","Magnetococcus_marinus_MC-1","Methylacidiphilum_infernorum_V4","Mycoplasma_genitalium_G37","Nostoc_punctiforme_PCC_73102","Opitutus_terrae_PB90-1","Pedobacter_heparinus_DSM_2366","Pirellula_staleyi_DSM_6068","Prochlorococcus_marinus_str._AS9601","Psychrobacter_arcticus_273-4","Rhizobium_leguminosarum_bv._trifolii_WSM1325_plasmid_pR132502","Rhodopirellula_baltica_SH_1_complete_genome","Rhodospirillum_rubrum_ATCC_11170","Rickettsia_rickettsii_str._Iowa","Shewanella_putrefaciens_CN-32","Synechococcus_elongatus_PCC_6301_DNA","Thermanaerovibrio_acidaminovorans_DSM_6589","Thermoanaerobacter_tengcongensis_MB4","Thermobaculum_terrenum_ATCC_BAA-798_chromosome_1","Thermodesulfovibrio_yellowstonii_DSM_11347","Thermomicrobium_roseum_DSM_5159_plasmid","Thermotoga_maritima_MSB8","Thermus_thermophilus_HB8_plasmid_pTT27_DNA"]
		names_archeas = ["Aeropyrum_pernix_K1_DNA","Archaeoglobus_fulgidus_DSM_4304","Archaeoglobus_profundus_DSM_5631","Caldivirga_maquilingensis_IC-167","Candidatus_Korarchaeum_cryptofilum_OPF8","Desulfurococcus_kamchatkensis_1221n","Haloarcula_marismortui_ATCC_43049_plasmid_pNG100","Halobacterium_sp._NRC-1","Halomicrobium_mukohataei_DSM_12286","Haloquadratum_walsbyi_DSM_16790_complete_genome","Halorhabdus_utahensis_DSM_12940","Halorubrum_lacusprofundi_ATCC_49239_chromosome_2","Haloterrigena_turkmenica_DSM_5511_plasmid_pHTUR04","Hyperthermus_butylicus_DSM_5456","Ignicoccus_hospitalis_KIN4_I","Metallosphaera_sedula_DSM_5348","Methanobrevibacter_ruminantium_M1","Methanobrevibacter_smithii_ATCC_35061","Methanocaldococcus_fervens_AG86_plasmid_pMEFER01","Methanocella_paludicola_SANAE_DNA","Methanococcoides_burtonii_DSM_6242","Methanococcus_aeolicus_Nankai-3","Methanococcus_maripaludis_C6","Methanococcus_vannielii_SB","Methanocorpusculum_labreanum_Z","Methanoculleus_marisnigri_JR1","Methanohalophilus_halophilus_strain_Z-7982_chromosome","Methanopyrus_kandleri_AV19","Methanoregula_boonei_6A8_chromosome","Methanosaeta_thermophila_PT","Methanosarcina_acetivorans_str._C2A","Methanosarcina_barkeri_str._Fusaro","Methanosarcina_mazei_strain_Goe1","Methanosphaera_stadtmanae_DSM_3091","Methanosphaerula_palustris_E1-9c_chromosome","Methanospirillum_hungatei_JF-1","Nanoarchaeum_equitans_Kin4-M","Natronomonas_pharaonis_DSM_2160_plasmid_PL131_complete_genome","Nitrosopumilus_maritimus_SCM1","Picrophilus_torridus_DSM_9790","Pyrobaculum_aerophilum_str._IM2","Pyrobaculum_arsenaticum_DSM_13514","Pyrobaculum_neutrophilum_V24Sta_chromosome","Pyrococcus_abyssi_complete_genome","Pyrococcus_furiosus_DSM_3638","Pyrococcus_horikoshii_OT3_DNA","Staphylothermus_marinus_F1","Sulfolobus_acidocaldarius_DSM_639","Sulfolobus_solfataricus_P2","Thermococcus_gammatolerans_EJ3","Thermofilum_pendens_Hrk_5","Thermoplasma_acidophilum_DSM_1728_complete_genome","Thermoplasma_volcanium_GSS1_DNA"]
		np.savetxt("distancematrix.csv", matrix, header = ",".join(names_archeas),delimiter = ",") #A transposer avec le code awk et a rajouter des headers (noms espaces separés par , et avec une , avant premier nom), pour pouvoir faire NJ
		
	