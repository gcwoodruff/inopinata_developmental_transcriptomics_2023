#count species-specific amino acid changes
import csv
import os, subprocess
import itertools

#set the directories

input_dir = '/home/gcw/Desktop/genome/species-specific_amino_acid_changes_12-22-17/08_transpose_folder_out/'

output_dir = '/home/gcw/Desktop/genome/species-specific_amino_acid_changes_12-22-17/09_species-specific_amino_acid_changes_b_out/'


#start the loop

for fil in os.listdir(input_dir):

	#open the file
	input_path = os.path.join(input_dir, fil)
	in_fil_open = open(input_dir + fil)
	in_fil_csv = csv.reader(in_fil_open)
		
	#start the species-specific residue lists

	brenneri_specific_residues = []
	briggsae_specific_residues = []
	doughertyi_specific_residues = []
	elegans_specific_residues = []
	inopinata_specific_residues = []
	latens_specific_residues = []
	nigoni_specific_residues = []
	remanei_specific_residues = []
	sinica_specific_residues = []
	sp26_specific_residues = []
	sp40_specific_residues = []
	tropicalis_specific_residues = []
	wallacei_specific_residues = []
	kamaaina_specific_residues = []
		

		#get the species-specific residues and add them to the lists

	for row in in_fil_csv:
		if row[0] != row[1] and row[1] == row[2] and row[1] == row[3] and row[1] == row[4] and row[1] == row[5] and row[1] == row[6] and row[1] == row[7] and row[1] == row[8] and row[1] == row[9] and row[1] == row[10] and row[1] == row[11] and row[1] == row[12] and row[1] == row[13]:
			brenneri_specific_residues.append(row)		
		if row[0] != row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			briggsae_specific_residues.append(row)		
		if row[0] == row[1] and row[0] != row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			doughertyi_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] != row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			elegans_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] != row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			inopinata_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] != row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			latens_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] != row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			nigoni_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] != row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			remanei_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] != row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			sinica_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] != row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			sp26_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] != row[10] and row[0] == row[11] and row[0] == row[12] and row[0] == row[13]:
			sp40_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] != row[11] and row[0] == row[12] and row[0] == row[13]:
			tropicalis_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] != row[12] and row[0] == row[13]:
			wallacei_specific_residues.append(row)		
		if row[0] == row[1] and row[0] == row[2] and row[0] == row[3] and row[0] == row[4] and row[0] == row[5] and row[0] == row[6] and row[0] == row[7] and row[0] == row[8] and row[0] == row[9] and row[0] == row[10] and row[0] == row[11] and row[0] == row[12] and row[0] != row[13]:
			kamaaina_specific_residues.append(row)		
	#count the number of species-specific residues
	number_of_brenneri_specific_residues = len(brenneri_specific_residues)
	number_of_briggsae_specific_residues = len(briggsae_specific_residues)
	number_of_doughertyi_specific_residues = len(doughertyi_specific_residues)
	number_of_elegans_specific_residues = len(elegans_specific_residues)
	number_of_inopinata_specific_residues = len(inopinata_specific_residues)
	number_of_latens_specific_residues = len(latens_specific_residues)
	number_of_nigoni_specific_residues = len(nigoni_specific_residues)
	number_of_remanei_specific_residues = len(remanei_specific_residues)
	number_of_sinica_specific_residues = len(sinica_specific_residues)
	number_of_sp26_specific_residues = len(sp26_specific_residues)
	number_of_sp40_specific_residues = len(sp40_specific_residues)
	number_of_tropicalis_specific_residues = len(tropicalis_specific_residues)
	number_of_wallacei_specific_residues = len(wallacei_specific_residues)
	number_of_kamaaina_specific_residues = len(kamaaina_specific_residues)

	#string-ify

	number_of_brenneri_specific_residues_str = str(number_of_brenneri_specific_residues)
	number_of_briggsae_specific_residues_str = str(number_of_briggsae_specific_residues)
	number_of_doughertyi_specific_residues_str = str(number_of_doughertyi_specific_residues)
	number_of_elegans_specific_residues_str = str(number_of_elegans_specific_residues)
	number_of_inopinata_specific_residues_str = str(number_of_inopinata_specific_residues)
	number_of_latens_specific_residues_str = str(number_of_latens_specific_residues)
	number_of_nigoni_specific_residues_str = str(number_of_nigoni_specific_residues)
	number_of_remanei_specific_residues_str = str(number_of_remanei_specific_residues)
	number_of_sinica_specific_residues_str = str(number_of_sinica_specific_residues)
	number_of_sp26_specific_residues_str = str(number_of_sp26_specific_residues)
	number_of_sp40_specific_residues_str = str(number_of_sp40_specific_residues)
	number_of_tropicalis_specific_residues_str = str(number_of_tropicalis_specific_residues)
	number_of_wallacei_specific_residues_str = str(number_of_wallacei_specific_residues)
	number_of_kamaaina_specific_residues_str = str(number_of_kamaaina_specific_residues)
	#open out file
	fil_out = output_dir + fil
	out_fil = open(fil_out, "w")
	#write it
	out_fil.write(fil + "\t" +
		number_of_brenneri_specific_residues_str + "\t" + number_of_briggsae_specific_residues_str + "\t" + number_of_doughertyi_specific_residues_str + "\t" + number_of_elegans_specific_residues_str + "\t" + number_of_inopinata_specific_residues_str + "\t" + number_of_latens_specific_residues_str + "\t" + number_of_nigoni_specific_residues_str + "\t" + 
		number_of_remanei_specific_residues_str + "\t" + number_of_sinica_specific_residues_str + "\t" + number_of_sp26_specific_residues_str + "\t" + number_of_sp40_specific_residues_str + "\t" + number_of_tropicalis_specific_residues_str + "\t" + number_of_wallacei_specific_residues_str + "\t" + number_of_kamaaina_specific_residues_str)
	out_fil.close()
		
		
		
		
		



