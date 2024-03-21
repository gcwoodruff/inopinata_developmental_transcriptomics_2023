import re
import os, subprocess
import itertools

input_dir = '/home/gcw/Desktop/genome/species-specific_amino_acid_changes_12-22-17/05_sed_4/'

output_dir = '/home/gcw/Desktop/genome/species-specific_amino_acid_changes_12-22-17/06_add_commas_folder_out/'

for fil in os.listdir(input_dir):
	fil_out = output_dir + fil


	with open(input_dir + fil, "r") as infile, open(fil_out, "w") as output_file:
	
		for line in infile:
			new_line = re.sub(r'([A-Z])', r'\1,', line)
			output_file.write(new_line)
	
	output_file.close()
