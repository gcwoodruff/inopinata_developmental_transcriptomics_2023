#getting transcripts of genes that align to transposons!

	#put working directory here
wkdir='/projects/phillipslab/gavincw/inopinata_RNAseq_1-20-20/'

mkdir 13_transposon_aligning_proteins_3-14-20

cd 13_transposon_aligning_proteins_3-14-20

mkdir 00_blast

#transposon db

	#all these files are from the paper Woodruff and Anastasia 2020 MBE

# /projects/phillipslab/gavincw/repeats_12-18-18/30_blastp/01_makeblastdb/transposon_db.pep

#elegans proteins

# $wkdir/03_worm_references/elegans/elegans_protein_longest_isoform.fa

#inopinata proteins

# $wkdir/03_worm_references/inopinata/inopinata_protein_longest_isoform.fa


#blast elegans proteins to transposon db

#module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 NCBI-Toolkit/18.0.0

blastp -db /projects/phillipslab/gavincw/repeats_12-18-18/30_blastp/01_makeblastdb/transposon_db -outfmt 6 -query elegans_protein_longest_isoform.fa -out $wkdir/13_transposon_aligning_proteins_3-14-20/00_blast/elegans.blastp_out -evalue 0.005

#blast inopinata proteins to transposon db


blastp -db /projects/phillipslab/gavincw/repeats_12-18-18/30_blastp/01_makeblastdb/transposon_db -outfmt 6 -query $wkdir/03_worm_references/inopinata/inopinata_protein_longest_isoform.fa -out $wkdir/13_transposon_aligning_proteins_3-14-20/00_blast/inopinata.blastp_out -evalue 0.005


mkdir $wkdir/13_transposon_aligning_proteins_3-14-20/01_cut_unique/

cd  $wkdir/13_transposon_aligning_proteins_3-14-20/00_blast/

#get unique hits
for i in *; do cut -f1 $i | uniq >  $wkdir/13_transposon_aligning_proteins_3-14-20/01_cut_unique/${i%.blast_out}; done


mkdir $wkdir/13_transposon_aligning_proteins_3-14-20/02_blast/



#get all the elegans protein ids
grep ">" $wkdir/03_worm_references/elegans/elegans_protein_longest_isoform.fa | sed -e 's/>//g' | sed -e 's/ wormpep.*//g' | sort | uniq > $wkdir/13_transposon_aligning_proteins_3-14-20/02_blast/00_tot_elegans_prot_id

#get unique elegans protein ids that align to transposons
sort $wkdir/13_transposon_aligning_proteins_3-14-20/01_cut_unique/elegans.blastp_out | uniq > $wkdir/13_transposon_aligning_proteins_3-14-20/02_blast/01_elegans_align_transposon_prot_id

#get the elegans proteins that DO NOT align to transposons

cd  $wkdir/13_transposon_aligning_proteins_3-14-20/02_blast/

comm -23 00_tot_elegans_prot_id 01_elegans_align_transposon_prot_id > 02_elegans_non-transposon_prot_id

#just get those sequences


perl fasta_filter.pl 02_elegans_non-transposon_prot_id elegans_protein_longest_isoform.fa 03_elegans_non-transposon_prot.fa 


sed -i -e 's/ wormpep.*//g' 03_elegans_non-transposon_prot.fa

sed -e '/^>/s/$/@/' -e 's/^>/#/' 03_elegans_non-transposon_prot.fa |\
tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
sort -u -t ' ' -f -k1,1 |\
sed -e 's/^/>/' -e 's/\t/\n/' > 04_elegans_non-transposon_prot.fa

#make a db with the elegans proteins that _do not_ align to transposon domains
makeblastdb -in  04_elegans_non-transposon_prot.fa -parse_seqids -out elegans_non-transposon_prot_db -title elegans_non-transposon_prot_db -dbtype prot


#get the inopinata transposon-aligning proteins
perl fasta_filter.pl $wkdir/13_transposon_aligning_proteins_3-14-20/01_cut_unique/inopinata.blastp_out $wkdir/03_worm_references/inopinata/inopinata_protein_longest_isoform.fa 05_sp34_transposon_prot.fa


#which elegans non-transposon-aligning genes align to the inopinata transposon-aligning genes?
blastp -db elegans_non-transposon_prot_db -outfmt 6 -query 05_sp34_transposon_prot.fa -out 06_sp34_transposon_prot_to_elegans_blast_out -evalue 0.005

#get the ids of inopinata-transposon-aligning genes that align to elegans-non-transposons

cut -f1 06_sp34_transposon_prot_to_elegans_blast_out | sort | uniq > 07_34_transposon_prot_align_to_elegans_list

#extract only those inopinata transposon-aligning genes that only align to transposons and not to other non-transposon elegans sequences!
comm -23 <(grep -Po '\S+' $wkdir/13_transposon_aligning_proteins_3-14-20/01_cut_unique/inopinata.blastp_out | sort) <(grep -Po '\S+' 07_34_transposon_prot_align_to_elegans_list | sort)  > 08_inopinata_transposon_protein-coding_genes

	#that's the list

grep -f 07_34_transposon_prot_align_to_elegans_list $wkdir/12_transposon_aligning_proteins/rlog_transformed_transcript_counts.tsv | 
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"transposon-aligning"}' > rlog_transformed_transcript_counts_TRANSPOSON_ALIGNING


grep -v -f 07_34_transposon_prot_align_to_elegans_list $wkdir/12_transposon_aligning_proteins/rlog_transformed_transcript_counts.tsv | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"NOT transposon-aligning"}' > rlog_transformed_transcript_counts_NOT_TRANSPOSON_ALIGNING
	#this was used for analysis and figures
cat rlog_transformed_transcript_counts_NOT_TRANSPOSON_ALIGNING rlog_transformed_transcript_counts_TRANSPOSON_ALIGNING > rlog_transformed_transcript_counts_transposon_tags.tsv
