#Getting species-specific amino acid replacements among Caenorhabditis species.

#For Woodruff et al. "Widespread changes in gene expression accompany body size evolution in nematodes"

#This work was done a while ago, in 2017. I am doing my best to cobble together the pieces of this workflow.

#These protein sets were retrieved:

#caenorhabditis_brenneri.PRJNA20035.WBPS9.protein.fa
#caenorhabditis_briggsae.PRJNA10731.WBPS9.protein.fa
#caenorhabditis_elegans.PRJNA13758.WBPS9.protein.fa
#Caenorhabditis_remanei_PX439_v1.proteins.fa
#caenorhabditis_sinica.PRJNA194557.WBPS9.protein.fa
#caenorhabditis_tropicalis.PRJNA53597.WBPS9.protein.fa
#CDOUG.caenorhabditis_doughertyi.JU1771.CGP2.proteins.fa
#Caenorhabditis_sp26_JU2190_v1.proteins.fa
#CSP40.caenorhabditis_sp40.JU2818.CGP2.proteins.fa
#wallacei_gen.vs.wallacei_cDNA_2016.11.02.complete.pep.fa
#nigoni.pc_gen.vs.nigoni_cDNA_2016.11.02.complete.pep.fa
#Caenorhabditis_latens_PX534_v1.proteins.fa
#Caenorhabditis_sp34_NK74SC_v710.proteins.fa
#and C. kamaaina

#New versions of these protein sets can be retrieved from WormBase ParaSite and/or caenorhabditis.org

#Largest isoforms were extracted and these protein sets were then used in an all-by-all blastp, files prepared by orthofinder

wkdir = 'path/to/working/directory'

cd $wkdir/orthofinder_prepare_blast_files/

./orthofinder -f $wkdir/caenorhabditis_proteins_longest_isoform/ -op

#blastp was done with these files using these commands: blastp_commands.sh , at https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/blastp_commands.sh

#version of blast: blast/2.2.30+

#run orthofinder on the blast output with only groups option

./orthofinder -b $wkdir/orthofinder_blast_output/ -og

#this generated a file, Orthogroups.csv (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/Orthogroups.csv)

#For earlier versions of Orthofinder, it did not generate output files with the number of proteins per species per orthogroup. To do this in 2017, I split the columns of the Orthogroups.csv into separate files-- each file having a column with the Orthogroup ID and a column of the proteins in that orthogroup for a given species. Then I used commas.py (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/commas.py) to count the commas in each row, which can be used to get the protein count in the orthogroup.

#This generated the file, orthogroup_counts.csv https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/orthogroup_counts.csv

#to get the single-copy orthogroups

awk 'BEGIN {FS="\t"} $2 == 1 && $3 == 1 && $4 == 1 && $5 == 1 && $6 == 1 && $7 == 1 && $8 == 1 && $9 == 1 && $10 == 1 && $11 == 1 && $12 == 1 && $13 == 1 && $14 == 1 {print} ' $wkdir/orthogroup_counts/orthogroup_counts.csv > $wkdir/orthogroup_counts/orthogroup_counts/1-to-1_orthogroups.tsv

#to get the single-copy orthogroup OG ID's

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' $wkdir/orthogroup_counts/orthogroup_counts/1-to-1_orthogroups.tsv > $wkdir/orthogroup_counts/orthogroup_counts/single-copy_orthogroup_OG_IDs.txt

#get the sequence id's of the genes in the single-copy orthogroups

grep -w -f single-copy_orthogroup_OG_IDs.txt Orthogroups.csv > single-copy_orthogroup_protein_ids.tsv

#then, fasta files containing the proteins for each single-copy orthologous group were generated. 
#This was ultimately done by transposing single-copy_orthogroup_protein_ids.tsv and separating it into 3,290 files listing the names of the proteins. These are in the file: prot_lists.zip (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/prot_lists.zip)

#Then, ALL of the longest-isoform Caenorhabditis proteins for the species of interest were put into a single file, and the protein sequences themselves were retrieved with fasta_filter.pl (Thanks Kevin Nyberg, who wrote this many years ago). This generated 3,290 fasta files, each containing 14 orthologous proteins, each from 14 species. These are here in the file fasta_filter_out.zip (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/fasta_filter_out.zip)

#These were then aligned with mafft

cd $wkdir/fasta_filter_out/

for i in *; do mafft --auto $i > $wkdir/mafft_out/$i; done

#These alignments were then trimmed with trimal

cd $wkdir/mafft_out/

for i in *.fa; do trimal -in $i -out $wkdir/trimal_out/$i -gt 1; done

	#these trimmed alignments are in trimal_out.zip ; 

#then, I apparently transformed these alignments into csv's all with an equal number of columns to ultimately count the number of species-specific amino acid replacements


cd $wkdir

mkdir 01_tr
mkdir 02_sed_1
mkdir 03_sed_2
mkdir 04_sed_3
mkdir 05_sed_4
mkdir 06_add_commas_folder_out
mkdir 07_sed_5
mkdir 08_transpose_folder_out
mkdir 09_species-specific_amino_acid_changes_b_out

cd $wkdir/trimal_out/

#replace all return with nothing
for i in *; do tr '\n' ' ' < $i >  $wkdir/01_tr/$i; done

cd $wkdir/01_tr/

#replace all > with \n>

for i in *; do sed -e 's/>/\n>/g' $i > $wkdir/02_sed_1/$i; done

cd $wkdir/02_sed_1

#remove all spaces

for i in *; do sed -e 's/\s//g' $i > $wkdir/03_sed_2/$i; done

cd $wkdir/03_sed_2
#remove first return
for i in *; do sed '1d' $i > $wkdir/04_sed_3/$i; done

cd $wkdir/04_sed_3
#remove non-AA stuff

for i in *; do sed -e 's/^.*bp//g' $i > $wkdir/05_sed_4/$i; done

cd $wkdir/05_sed_4

#add commas with add_commas_folder.py

cd ..

python add_commas_folder.py


#after add_commas_folder.py, do...

cd $wkdir/06_add_commas_folder_out

for i in *; do sed -e 's/,$//g' $i > $wkdir/07_sed_5/$i; done

#transpose alignment csv's

cd ..

python transpose_folder.py

#now, for the moment, species-specific_amino_acid_changes_b.py

python species-specific_amino_acid_changes_b.py

#this made 3,290 files, each being a row of the species-specific amino acid replacement data. These were then combined

cd $wkdir/species-specific_amino_acid_changes_b_py_out/

cat * > $wkdir/number of species-specific_aa_changes_for_single-copy_orthologs.tsv

#this file was opened in LibreOffice and exported as a csv: This is called number_of_species-specific_aa_changes_for_single_copy_orthologs_b.csv (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/number_of_species-specific_aa_changes_for_single_copy_orthologs_b.csv)

#This was used as data for further analysis, statistics, and visualization in species-specific_amino_acid_replacements_notes.R

#also, I made a key for C. elegans protein/gene ids

cd $wkdir/trimal_out/

mkdir $wkdir/eleg_gene_id/

#get the C. elegans id's. because trimal also prints the length of the sequence, this will also get the length of the alignment (aka prot_length)
for i in *; do grep ">CELEG" $i > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/AA_old_ubuntu_desktop/species-specific_amino_acid_changes_12-22-17/eleg_gene_id/$i; done

#combine all elegans gene id's
cd $wkdir/eleg_gene_id/

cat * > $wkdir/eleg_gene_id_plus_prot_length.txt

#get the list of the fasta file names
ls > $wkdir/fa_file_id.txt

cd $wkdir/
#combine C. elegans gene id's with fasta file id's to get the key
paste eleg_gene_id_plus_prot_length.txt fa_file_id.txt > eleg_gene_id_prot_length_fa_file_id.txt

#I used wormbase simplemine with this list of C. elegans gene id's to get the common/public names. This was added with LibreOffice. This is the file: eleg_gene_id_prot_length_fa_file_id_repaired.csv (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/eleg_gene_id_prot_length_fa_file_id_repaired.csv)

#Next, species-specific_amino_acid_replacements_notes.R

















