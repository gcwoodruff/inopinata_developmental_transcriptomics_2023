import csv
from itertools import izip
import os, subprocess

input_dir = '/home/gcw/Desktop/genome/species-specific_amino_acid_changes_12-22-17/07_sed_5/'

output_dir = '/home/gcw/Desktop/genome/species-specific_amino_acid_changes_12-22-17/08_transpose_folder_out/'



for fil in os.listdir(input_dir):

	fil_out = output_dir + fil

	#open the file
	input_path = os.path.join(input_dir, fil)
	in_fil_open = open(input_dir + fil)
	in_fil_csv = csv.reader(in_fil_open)


	a = izip(*in_fil_csv)
	csv.writer(open(fil_out, "wb")).writerows(a)

