#This is the code for generating INTERPRO domain annotations for Woodruff et al. "Widespread changes in gene expression accompany body size evolution in nematodes"

#installed interproscan-5.65-97.0-64-bit.tar.gz

#load the required module
#module load Java/11.0.16




#make some directories
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/briggsae
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/elegans
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/inopinata
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/nigoni
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/remanei

#run interproscan (protein sets retrieved/processed as described in pipeline.sh https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/pipeline.sh)
cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/interproscan/interproscan-5.65-97.0

./interproscan.sh -i /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/inopinata_RNA-seq_phillipslab/inopinata_RNAseq_1-20-20/11_new_ortho_fork/01_get_largest_isoform/elegans.fa -b /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/elegans/

./interproscan.sh -i /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/inopinata_RNA-seq_phillipslab/inopinata_RNAseq_1-20-20/11_new_ortho_fork/01_get_largest_isoform/inopinata.fa -b /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/inopinata/

./interproscan.sh -i /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/inopinata_RNA-seq_phillipslab/inopinata_RNAseq_1-20-20/11_new_ortho_fork/01_get_largest_isoform/briggsae.fa -b /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/briggsae/

./interproscan.sh -i /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/inopinata_RNA-seq_phillipslab/inopinata_RNAseq_1-20-20/11_new_ortho_fork/01_get_largest_isoform/nigoni.fa -b /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/nigoni/

./interproscan.sh -i /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/inopinata_RNA-seq_phillipslab/inopinata_RNAseq_1-20-20/11_new_ortho_fork/01_get_largest_isoform/remanei.fa -b /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024/00_interproscan/remanei/

#okay, cool interproscan.

#prep for R locally
	#get rows that have interpro id's and human-readable descriptions
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$12,$13}' inopinata.fa.tsv > inopinata_02.tsv
	#remove rows with no ids and descriptions
awk 'BEGIN {FS="\t"} $2 != "-" {print} ' inopinata_02.tsv > inopinata_03.tsv
	#remove duplicate entries -- NO DUPLICATE DOMAINS PER GENE
sort inopinata_03.tsv | uniq > inopinata_04.tsv


	#get rows that have interpro id's and human-readable descriptions
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$12,$13}' elegans.fa.tsv > elegans_02.tsv
	#remove rows with no ids and descriptions
awk 'BEGIN {FS="\t"} $2 != "-" {print} ' elegans_02.tsv > elegans_03.tsv
	#remove duplicate entries -- NO DUPLICATE DOMAINS PER GENE
sort elegans_03.tsv | uniq > elegans_04.tsv








