#Trying to get the number of entries per column, to figure out the expansion/contraction of ortho-copies in caenorhabditis. This time, trying to get at empty/null cells too.

import csv

bre = open("brenneri_copies")
bri = open("briggsae_copies")
dou = open("doughertyi_copies")
ele = open("elegans_copies")
ino = open("inopinata_copies")
lat = open("latens_copies")
nig = open("nigoni_copies")
rem = open("remanei_copies")
sin = open("sinica_copies")
sp26 = open("sp26_copies")
sp40 = open("sp40_copies")
tro = open("tropicalis_copies")
wal = open("wallacei_copies")


bre_csv = csv.reader(bre, delimiter='\t')
bri_csv = csv.reader(bri, delimiter='\t')
dou_csv = csv.reader(dou, delimiter='\t')
ele_csv = csv.reader(ele, delimiter='\t')
ino_csv = csv.reader(ino, delimiter='\t')
lat_csv = csv.reader(lat, delimiter='\t')
nig_csv = csv.reader(nig, delimiter='\t')
rem_csv = csv.reader(rem, delimiter='\t')
sin_csv = csv.reader(sin, delimiter='\t')
sp26_csv = csv.reader(sp26, delimiter='\t')
sp40_csv = csv.reader(sp40, delimiter='\t')
tro_csv = csv.reader(tro, delimiter='\t')
wal_csv = csv.reader(wal, delimiter='\t')

for row in bre_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_bre = open("brenneri_gene_counts", "a")
	output_bre.write(row[0] + "\t" + gene_count_str + "\n")
	output_bre.close()

for row in bri_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_bri = open("briggsae_gene_counts", "a")
	output_bri.write(row[0] + "\t" + gene_count_str + "\n")
	output_bri.close()

for row in dou_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_dou = open("doughertyi_gene_counts", "a")
	output_dou.write(row[0] + "\t" + gene_count_str + "\n")
	output_dou.close()

for row in ele_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_ele = open("elegans_gene_counts", "a")
	output_ele.write(row[0] + "\t" + gene_count_str + "\n")
	output_ele.close()

for row in ino_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_ino = open("inopinata_gene_counts", "a")
	output_ino.write(row[0] + "\t" + gene_count_str + "\n")
	output_ino.close()

for row in lat_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_lat = open("latens_gene_counts", "a")
	output_lat.write(row[0] + "\t" + gene_count_str + "\n")
	output_lat.close()

for row in nig_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_nig = open("nigoni_gene_counts", "a")
	output_nig.write(row[0] + "\t" + gene_count_str + "\n")
	output_nig.close()

for row in rem_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_rem = open("remanei_gene_counts", "a")
	output_rem.write(row[0] + "\t" + gene_count_str + "\n")
	output_rem.close()

for row in sin_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_sin = open("sinica_gene_counts", "a")
	output_sin.write(row[0] + "\t" + gene_count_str + "\n")
	output_sin.close()

for row in sp26_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_sp26 = open("sp26_gene_counts", "a")
	output_sp26.write(row[0] + "\t" + gene_count_str + "\n")
	output_sp26.close()

for row in sp40_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_sp40 = open("sp40_gene_counts", "a")
	output_sp40.write(row[0] + "\t" + gene_count_str + "\n")
	output_sp40.close()

for row in tro_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_tro = open("tropicalis_gene_counts", "a")
	output_tro.write(row[0] + "\t" + gene_count_str + "\n")
	output_tro.close()

for row in wal_csv:
	genes = row[1]
	comma = genes.count(',')
	if genes in (None, ""):
		null = 1
	else: 
		null = 0
	gene_count = (comma + 1) - null
	gene_count_str = str(gene_count)
	output_wal = open("wallacei_gene_counts", "a")
	output_wal.write(row[0] + "\t" + gene_count_str + "\n")
	output_wal.close()

