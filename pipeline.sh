
#C. inopinata developmental transcriptomics, Woodruff et al. 2023
#gcwoodruff@ou.edu

#put working directory here
wkdir="/path/to/files/"

mkdir $wkdir/00_raw_read_links
#link raw reads here
cd $wkdir/00_raw_read_links

#FastQC was used to evaluate read quality.

mkdir $wkdir/01_FastQC/

cd $wkdir/00_raw_read_links

fastqc -o $wkdir/01_FastQC/ -f fastq Undetermined_S0_L001_I1_001.fastq.gz Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz
fastqc -o $wkdir/01_FastQC/ -f fastq Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz
	#the reads looked alright!

#demultiplex reads with stacks process_shortreads

mkdir $wkdir/02_stacks_process_short_reads
cd $wkdir/00_raw_read_links

process_shortreads -1 Undetermined_S0_L001_R1_001.fastq.gz -2 Undetermined_S0_L001_R2_001.fastq.gz -i gzfastq -b $wkdir/barcodes_2.txt -o $wkdir/02_stacks_process_short_reads/ -q -c -r --index_null



#get and prep worm references
mkdir $wkdir/03_worm_references
mkdir $wkdir/03_worm_references/elegans
mkdir $wkdir/03_worm_references/inopinata

cd $wkdir/03_worm_references/elegans

wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.genomic_softmasked.fa.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.mRNA_transcripts.fa.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.annotations.gff3.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.canonical_geneset.gtf.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.protein.fa.gz

gunzip *

cd $wkdir/03_worm_references/inopinata


wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_inopinata/PRJDB5687/c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_inopinata/PRJDB5687/c_inopinata.PRJDB5687.WS275.genomic_softmasked.fa.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_inopinata/PRJDB5687/c_inopinata.PRJDB5687.WS275.annotations.gff3.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_inopinata/PRJDB5687/c_inopinata.PRJDB5687.WS275.canonical_geneset.gtf.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_inopinata/PRJDB5687/c_inopinata.PRJDB5687.WS275.protein.fa.gz

gunzip *

#get longest isoform of each transcript

cd $wkdir/03_worm_references/elegans/

perl $wkdir/scripts/get_largest_isoforms.pl -i c_elegans.PRJNA13758.WS275.mRNA_transcripts.fa -t wb > elegans_mRNA_transcripts_longest_isoform.fa 

grep -c '^>' elegans_mRNA_transcripts_longest_isoform.fa
# 20040

perl $wkdir/scripts/get_largest_isoforms.pl -i c_elegans.PRJNA13758.WS275.protein.fa -t wb > elegans_protein_longest_isoform.fa 

grep -c '^>' elegans_protein_longest_isoform.fa
# 20040


cd  $wkdir/03_worm_references/inopinata

perl $wkdir/scripts/get_largest_isoforms.pl -i c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa -t aug_like > inopinata_mRNA_transcripts_longest_isoform.fa 


grep pseudo c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa

# >Sp34_10027510.pseudo1 gene=Sp34_10027510
# >Sp34_10090810.pseudo1 gene=Sp34_10090810
# >Sp34_30332210.pseudo1 gene=Sp34_30332210
# >Sp34_50073910.pseudo1 gene=Sp34_50073910
# >Sp34_50280940.pseudo1 gene=Sp34_50280940

#put this in inopinata_pseudogenes_list.txt

#remove pseudogenes
faSomeRecords -exclude c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa inopinata_pseudogenes_list.txt inopinata_mRNA_transcripts_no_pseudogenes.fa

#get longest isoform
perl $wkdir/scripts/get_largest_isoforms.pl -i inopinata_mRNA_transcripts_no_pseudogenes.fa -t aug_like > inopinata_mRNA_transcripts_longest_isoform.fa 

faSomeRecords -exclude c_inopinata.PRJDB5687.WS275.protein.fa inopinata_pseudogenes_list.txt inopinata_protein_no_pseudogenes.fa


perl $wkdir/scripts/get_largest_isoforms.pl -i inopinata_protein_no_pseudogenes.fa -t aug_like > inopinata_protein_longest_isoform.fa &

grep -c '^>' inopinata_mRNA_transcripts_longest_isoform.fa
# 21604
grep -c '^>' inopinata_protein_longest_isoform.fa
# 21442

grep "^>" inopinata_protein_longest_isoform.fa | sed -e 's/^>//g' > inopinata_protein_longest_isoform_list.txt


grep -c "^>" inopinata_mRNA_transcripts_no_pseudogenes.fa
# 21604
grep -c "^>" inopinata_protein_no_pseudogenes.fa
# 21442

grep -c "^>" c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa
#21609

grep -c "^>" c_inopinata.PRJDB5687.WS275.protein.fa
# 21443

grep -c "t1" c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa
# 21604

grep -c "t2" c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa
# 0

grep -c "pseudo" c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa
#5 ; no alternative transcripts here. why are so many proteins missing?

#translate the mrna to protein with EMBOSS

#module load easybuild  ifort/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 EMBOSS/6.6.0

transeq -sequence c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa -outseq c_inopinata.PRJDB5687.WS275.mRNA_translated_to_protein.fa

#will just use that as it doesn't seem like there are alternative transcripts here.


#make transcript indices
cd $wkdir/03_worm_references/elegans

mkdir salmon_index
cd salmon_index

salmon index -t $wkdir/03_worm_references/elegans/elegans_mRNA_transcripts_longest_isoform.fa -i elegans_mRNA_lngst_index 

cd $wkdir/03_worm_references/inopinata
mkdir salmon_index
cd salmon_index

salmon index -t $wkdir/03_worm_references/inopinata/c_inopinata.PRJDB5687.WS275.mRNA_transcripts.fa -i inopinata_mRNA_index 

#map reads/quantify transcript abundance with salmon!

#samples 1-15 are C. inopinata; samples 16-30 are C. elegans

mkdir $wkdir/04_salmon_quant

salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_2_all.1.fq.gz  -2 sample_2_all.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_2_all
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_3_all.1.fq.gz  -2 sample_3_all.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_3_all
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_4.1.fq.gz  -2 sample_4.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_4
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_5.1.fq.gz  -2 sample_5.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_5
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_6.1.fq.gz  -2 sample_6.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_6
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_7.1.fq.gz  -2 sample_7.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_7
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_8.1.fq.gz  -2 sample_8.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_8
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_9.1.fq.gz  -2 sample_9.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_9
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_10.1.fq.gz  -2 sample_10.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_10
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_11.1.fq.gz  -2 sample_11.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_11
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_12.1.fq.gz  -2 sample_12.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_12
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_13.1.fq.gz  -2 sample_13.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_13
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_14.1.fq.gz  -2 sample_14.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_14
salmon quant -i $wkdir/03_worm_references/inopinata/salmon_index/inopinata_mRNA_index -l A -1 sample_15.1.fq.gz  -2 sample_15.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_15
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_16.1.fq.gz  -2 sample_16.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_16
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_17.1.fq.gz  -2 sample_17.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_17
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_18.1.fq.gz  -2 sample_18.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_18
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_19.1.fq.gz  -2 sample_19.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_19
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_20.1.fq.gz  -2 sample_20.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_20
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_21.1.fq.gz  -2 sample_21.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_21
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_22.1.fq.gz  -2 sample_22.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_22
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_23.1.fq.gz  -2 sample_23.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_23
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_24.1.fq.gz  -2 sample_24.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_24
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_25.1.fq.gz  -2 sample_25.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_25
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_26.1.fq.gz  -2 sample_26.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_26
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_27.1.fq.gz  -2 sample_27.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_27
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_28.1.fq.gz  -2 sample_28.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_28
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_29.1.fq.gz  -2 sample_29.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_29
salmon quant -i $wkdir/03_worm_references/elegans/salmon_index/elegans_mRNA_lngst_index  -l A -1 sample_30.1.fq.gz  -2 sample_30.2.fq.gz -p 8 --validateMappings --gcBias -o $wkdir/04_salmon_quant/sample_30



#ok, get percent mapped for all reads

mkdir $wkdir/04_salmon_quant/percent_mapped
mkdir $wkdir/04_salmon_quant/percent_mapped/00_grep
mkdir $wkdir/04_salmon_quant/percent_mapped/01_awk

while read F; do
  grep 'percent_mapped' $wkdir/04_salmon_quant/$F/aux_info/meta_info.json > $wkdir/04_salmon_quant/percent_mapped/00_grep/$F
  cd $wkdir/04_salmon_quant/percent_mapped/00_grep/
  sed -i -e 's/.*: //g' $F
  sed -i -e 's/,//g' $F
  awk 'BEGIN {FS="\t"} {OFS="\t"} {print FILENAME,$0}' $F > $wkdir/04_salmon_quant/percent_mapped/01_awk/$F
done < $wkdir/sample_list.txt

	#huh, why does process shortreads log have retained reads for sample_8 but there are none in the fq.gz????
	#I re-did stacks process shortreads in folder 05_re-do_stacks_process_shortreads and got the same result, hmm.

	#anyway moving on for now.........

cd $wkdir/04_salmon_quant/percent_mapped/01_awk/
cat * > $wkdir/04_salmon_quant/percent_mapped/percent_mapped.txt

cat $wkdir/04_salmon_quant/percent_mapped/percent_mapped.txt

#sample_1	78.4897617091672
#sample_10	79.3086614340584
#sample_11	76.44558531861154
#sample_12	81.37446715719112
#sample_13	71.32483350589215
#sample_14	78.79329775791385
#sample_15	79.18796713352701
#sample_16	7.494076579282733
#sample_17	87.37381184009344
#sample_18	89.71910790402774
#sample_19	85.97463759842624
#sample_20	90.1255736661657
#sample_21	87.35546924900942
#sample_22	91.209596064287
#sample_23	91.56963599130725
#sample_24	90.96682312643665
#sample_25	86.4074215807347
#sample_26	91.56004172434503
#sample_27	92.56812361579876
#sample_28	93.16100328321116
#sample_29	92.99990352444957
#sample_2_all	77.20021513970106
#sample_30	93.02610764962603
#sample_3_all	79.78622801260113
#sample_4	78.70112424374432
#sample_5	75.03906200348989
#sample_6	77.54326382685348
#sample_7	74.85824004426374
#sample_9	76.5424546360237

#as sample_8 and sample_16 revealed very low number of mapped reads, they will be excluded from subsequent analyses





#getting 1-1 orthologs
mkdir $wkdir/11_new_ortho_fork

cd $wkdir/11_new_ortho_fork

mkdir 00_links
#elegans and inopinata proteins from the available genome assemblies
ln -s $wkdir/03_worm_references/elegans/c_elegans.PRJNA13758.WS275.protein.fa $wkdir/11_new_ortho_fork/00_links/elegans_prot.fa

ln -s $wkdir/03_worm_references/inopinata/c_inopinata.PRJDB5687.WS275.protein.fa $wkdir/11_new_ortho_fork/00_links/inopinata_prot.fa

#briggsae proteins
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/caenorhabditis_briggsae/PRJNA10731/caenorhabditis_briggsae.PRJNA10731.WBPS14.protein.fa.gz
#nigoni proteins
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/caenorhabditis_nigoni/PRJNA384657/caenorhabditis_nigoni.PRJNA384657.WBPS14.protein.fa.gz
#remanei proteins
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/183/535/GCA_010183535.1_CRPX506/GCA_010183535.1_CRPX506_protein.faa.gz

gunzip *.gz

mv caenorhabditis_briggsae.PRJNA10731.WBPS14.protein.fa briggsae_prot.fa
mv caenorhabditis_nigoni.PRJNA384657.WBPS14.protein.fa nigoni_prot.fa
mv GCA_010183535.1_CRPX506_protein.faa remanei_prot.fa

#it seems like remanei, inopinata do not have alternative transcripts
#elegans, briggsae, nigoni appear to have them...

#getting only one isoform per protein for orthofinder
cd ..
mkdir 01_get_largest_isoform

cd $wkdir/11_new_ortho_fork/00_links

perl $wkdir/scripts/get_largest_isoforms.pl -i elegans_prot.fa -t wb > $wkdir/11_new_ortho_fork/01_get_largest_isoform/elegans_prot.fa &

perl $wkdir/scripts/get_largest_isoforms.pl -i briggsae_prot.fa -t wb > $wkdir/11_new_ortho_fork/01_get_largest_isoform/briggsae_prot.fa &

perl $wkdir/scripts/get_largest_isoforms.pl -i nigoni_prot.fa -t col3 > $wkdir/11_new_ortho_fork/01_get_largest_isoform/nigoni_prot.fa &


#

cp $wkdir/03_worm_references/inopinata/inopinata_protein_longest_isoform.fa  $wkdir/11_new_ortho_fork/01_get_largest_isoform/inopinata_prot.fa

cp $wkdir/11_new_ortho_fork/00_links/remanei_prot.fa $wkdir/11_new_ortho_fork/01_get_largest_isoform/remanei_prot.fa

# [gavincw@n120 01_get_largest_isoform]$ for i in *; do grep -c "^>" $i; done
# 20829
# 20040
# 21442
# 29167
# 26189

cd $wkdir/11_new_ortho_fork/01_get_largest_isoform/

mv briggsae_prot.fa briggsae.fa
mv elegans_prot.fa elegans.fa
mv inopinata_prot.fa inopinata.fa
mv nigoni_prot.fa nigoni.fa
mv remanei_prot.fa remanei.fa

#blast for orthofinder
mkdir $wkdir/11_new_ortho_fork/02_orthofinder_prep_blast

cd $wkdir/11_new_ortho_fork/02_orthofinder_prep_blast

#module load easybuild BLAST+/2.7.1

orthofinder -op -S blast -f $wkdir/11_new_ortho_fork/01_get_largest_isoform/ &

#here are the prepped blast commands

blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species3.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies3 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast3_3.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species4.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies3 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast4_3.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species3.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies4 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast3_4.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species4.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies4 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast4_4.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species3.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies2 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast3_2.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species2.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies3 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast2_3.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species3.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies0 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast3_0.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species0.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies3 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast0_3.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species3.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies1 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast3_1.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species1.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies3 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast1_3.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species4.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies2 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast4_2.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species2.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies4 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast2_4.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species4.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies0 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast4_0.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species0.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies4 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast0_4.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species4.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies1 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast4_1.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species1.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies4 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast1_4.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species2.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies2 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast2_2.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species2.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies0 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast2_0.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species0.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies2 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast0_2.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species0.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies0 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast0_0.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species2.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies1 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast2_1.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species1.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies2 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast1_2.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species1.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies0 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast1_0.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species0.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies1 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast0_1.txt
blastp -outfmt 6 -evalue 0.001 -query $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Species1.fa -db $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/BlastDBSpecies1 -out $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/Blast1_1.txt


#then run orthofinder

orthofinder -S blast -M msa -b $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory -a 10

#####okay, do some of this again

mkdir $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs

#get the elegans-inopinata 1-1 orthologs
awk 'BEGIN {FS="\t"} {OFS="\t"} $3 == 1 && $4 == 1 {print} ' $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/OrthoFinder/Results_Feb20/Orthogroups/Orthogroups.GeneCount.tsv | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$3,$4}'  > $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/00_1-1_orthogroups

#get just the OrthoFinder OG id's
cd $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}'  00_1-1_orthogroups > 01_1-1_orthogroups
#get elegans and inopinata ID's
LC_ALL=C fgrep -w -f 01_1-1_orthogroups $wkdir/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/OrthoFinder/Results_Feb20/Orthogroups/Orthogroups.tsv | awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$3,$4}' > 02_1-1_orthogroups

wc -l 02_1-1_orthogroups

# 10719 02_1-1_orthogroups woooo a bunch more!


#but, these are the old-school elegans sequence id's.. I want human-readable ID's! Now, we go on a journey...
#just get the elegans genes with 1-1 inopinata orthologs
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2}' 02_1-1_orthogroups > 03_elegans_gene_ids

#remove alternative transcript identifiers
sed -i -e 's/[a-z]$//g' 03_elegans_gene_ids

#
#I did a simplemine search for all elegans genes earlier
	#see this tool to get simplemine_results.txt ()
cp $wkdir/06_get_1-1_inopinata_elegans_orthologs/simplemine_results.txt simplemine_results.txt

#human-redable gene names are here
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5}' simplemine_results.txt > 04_simplemine_four_cols

#let's connect with the elegans genes with orthologs!
LC_ALL=C fgrep -w -f 03_elegans_gene_ids 04_simplemine_four_cols > 05_1-1_og_simplemine_ids &


wc -l 05_1-1_og_simplemine_ids

# 10714 05_1-1_og_simplemine_ids
	#um, alright we are only missing four!!!
	#what's going on???

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4}' 05_1-1_og_simplemine_ids > 06_elg_seq_id

	#these are my genes with no simplemine entries
comm -23 <(grep -Po '\S+' 03_elegans_gene_ids | sort) <(grep -Po '\S+' 06_elg_seq_id | sort) > 07_no_simplemine_og

grep D1081.20 05_1-1_og_simplemine_ids

#found these on wormbase, will add

#WBGene00304985	D1081.20	Live	D1081.20
#WBGene00304992	F35C12.7	Live	F35C12.7
#WBGene00304816	R11G1.11	Live	R11G1.11

#cat 

cat 04_simplemine_four_cols 08_extra_simplemine_lines_2-21-20 > 09_simplemine_four_cols

#doing it again with revised simplemine file
LC_ALL=C fgrep -w -f 03_elegans_gene_ids 09_simplemine_four_cols > 10_1-1_og_simplemine_ids &


wc -l 10_1-1_og_simplemine_ids

#10717 ; that's fine

#elegans ids
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $4,$2}' 10_1-1_og_simplemine_ids > 11_elg_seq_ids
#elegans id and inopinata id
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3}' 02_1-1_orthogroups > 12_elg_ino_ortho_seq_ids


#remove alternative transcript identifiers

sed -i -e 's/[a-z]\t/\t/g' 12_elg_ino_ortho_seq_ids

#merge

######
#This to line 306 is just a matter of getting the 1-1 ortholog gene ID file merged 

perl /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl -k -e "no_key" 11_elg_seq_ids 12_elg_ino_ortho_seq_ids  2> merge_pl.error > 13_merge_og

cat cat merge_pl.error

# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'Y105E8A.7' already exists for file '11_elg_seq_ids'. Not adding value 'lev-10'.
# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'ZC416.8' already exists for file '11_elg_seq_ids'. Not adding value 'unc-17'.
# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'B0564.1' already exists for file '11_elg_seq_ids'. Not adding value 'tin-9.2'.
# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'ZC416.8' already exists for file '12_elg_ino_ortho_seq_ids'. Not adding value 'Sp34_40085300.t1'.

grep "no_key" 13_merge_og

# F54H12.10	no_key	Sp34_30203000.t1
# T24H7.9	no_key	Sp34_20128400.t1
# Y80D4G.2	no_key	Sp34_30047220.t1
# B0252.11	no_key	Sp34_20266820.t1

#first figure out no key things

cd $wkdir/11_new_ortho_fork/00_links

grep F54H12.10 elegans_prot.fa

# >F54H12.10 wormpep=CE53828 gene=WBGene00305131 status=Confirmed


cd $wkdir/11_new_ortho_fork/00_links

grep T24H7.9 elegans_prot.fa

#>T24H7.9 wormpep=CE21205 gene=WBGene00305157 status=Confirmed uniprot=Q7JPE4 insdc=CCD66148.1 product="GyrI-like domain-containing protein"


grep Y80D4G.2 elegans_prot.fa

#>Y80D4G.2 wormpep=CE53928 gene=WBGene00305069 status=Confirmed insdc=CAA005918.1

grep B0252.11 elegans_prot.fa
#>B0252.11 wormpep=CE53867 gene=WBGene00305086 status=Confirmed



cd $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs

grep WBGene00305131 simplemine_results.txt
#nothing.....
#weird

grep CE53828 simplemine_results.txt
#nothing  also this gene is not on wormbase... how????

#ok moving on

grep WBGene00305157 simplemine_results.txt
grep WBGene00305069 simplemine_results.txt
grep WBGene00305086 simplemine_results.txt
#all nothing

grep CE21205 simplemine_results.txt
grep CE53928 simplemine_results.txt
grep CE53867 simplemine_results.txt

#ok, yep not in simplemine and also not in WB, so weird. will manually edit these.


sed -i -e 's/F54H12.10\tno_key/F54H12.10\tF54H12.10/g' 13_merge_og
sed -i -e 's/T24H7.9\tno_key/T24H7.9\tT24H7.9/g' 13_merge_og
sed -i -e 's/Y80D4G.2\tno_key/Y80D4G.2\tY80D4G.2/g' 13_merge_og
sed -i -e 's/B0252.11\tno_key/B0252.11\tB0252.11/g' 13_merge_og

#ok, now the duplicates...



# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'Y105E8A.7' already exists for file '11_elg_seq_ids'. Not adding value 'lev-10'.
# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'ZC416.8' already exists for file '11_elg_seq_ids'. Not adding value 'unc-17'.
# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'B0564.1' already exists for file '11_elg_seq_ids'. Not adding value 'tin-9.2'.
# /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl: Key 'ZC416.8' already exists for file '12_elg_ino_ortho_seq_ids'. Not adding value 'Sp34_40085300.t1'.


grep Y105E8A.7 13_merge_og
grep ZC416.8 13_merge_og
grep B0564.1 13_merge_og
grep ZC416.8 13_merge_og


# grep Y105E8A.7 11_elg_seq_ids
# grep ZC416.8 11_elg_seq_ids
# grep B0564.1 11_elg_seq_ids
# grep ZC416.8 11_elg_seq_ids
# ok, these do have duplicate entries in WB which is annoying ; will just take those for now. if these genes arise again... take notice
	#and, right, at least unc-17/cha-1 is a nested gene. so for legacy purposes this makes sense.


grep ZC416.8 12_elg_ino_ortho_seq_ids

#ZC416.8	Sp34_40085400.t1
#ZC416.8	Sp34_40085300.t1


#this seems to be a bigger error, ok


grep ZC416.8  02_1-1_orthogroups
	#ok, we messed up here with respect to the longest isoform stuff? somehow that script left this there. ok, let's look at the final merge

grep ZC416.8 13_merge_og

# ZC416.8	cha-1	Sp34_40085400.t1

#this is fine

#ok, let's extract elegans and inopinata seq ID's

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' 13_merge_og > 14_elegans_seq_id_og
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $3}' 13_merge_og > 15_inopinata_seq_id_og



#connect previous salmon quants with the orthologs
mkdir $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs
cd 03_salmon_quant_1-1_orthologs

mkdir 00_links
cd 00_links
mkdir elegans
mkdir inopinata



#elegans samples links

cd $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/

# LC_ALL=C fgrep -w -f 14_elegans_seq_id_og 04_simplemine_four_cols > 16_1-1_og_simplemine_ids &

while read F; do
  ln -s $wkdir/04_salmon_quant/$F/quant.sf $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/00_links/elegans/$F
done < $wkdir/elegans_samples.txt


#inopinata samples links

while read F; do
  ln -s $wkdir/04_salmon_quant/$F/quant.sf  $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/00_links/inopinata/$F
done < $wkdir/inopinata_samples.txt

#remove bad samples

cd inopinata
rm sample_8

cd ..
cd elegans
rm sample_16

#
cd ..
cd ..
mkdir 01_grep
cd 01_grep
mkdir inopinata
mkdir elegans

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/00_links/elegans/


#remove ".1"

for i in *; do perl -p -i -e "s/\.[0-9]+\t/\t/" $i; done &


#remove lowercase letter

for i in *; do sed -i -e 's/[a-z]\t/\t/1' $i; done &


#get 1-1 orthologs

for i in *; do LC_ALL=C fgrep -w -f $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/14_elegans_seq_id_og $i > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/elegans/$i; done &


cd  $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/elegans
for i in *; do wc -l $i; done
	#10721

wc -l $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/14_elegans_seq_id_og
	#10718
	#ok... it's close let's just keep going!!!





cd  $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/00_links/inopinata/

for i in *; do LC_ALL=C fgrep -w -f  $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/15_inopinata_seq_id_og $i > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/inopinata/$i; done &

cd  $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/inopinata
for i in *; do wc -l $i; done
	#10718

wc -l $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/15_inopinata_seq_id_og
	#10718 , sweet, just use merge to remove the things that are not shared?

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/elegans/
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' sample_17 > 00_elg_trans_id


comm -23 <(grep -Po '\S+'  00_elg_trans_id | sort) <(grep -Po '\S+' $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/14_elegans_seq_id_og | sort) > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/extra_elg_transc

grep -f $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/extra_elg_transc $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/inopinata/sample_1
	#cool, just _remove_ those lines!

mkdir $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/02_grep_remove_extra_transc_elegans/

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/elegans/
rm 00_elg_trans_id
for i in *; do LC_ALL=C fgrep -w -v -f  $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/extra_elg_transc $i > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/02_grep_remove_extra_transc_elegans/$i; done &

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/02_grep_remove_extra_transc_elegans/

for i in *; do wc -l $i; done



grep -f $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/extra_elg_transc sample_1
	#these are duplicate entries, okay.

grep -f $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/extra_elg_transc $wkdir/04_salmon_quant/sample_17/quant.sf

grep B0564.11 $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/02_grep_remove_extra_transc_elegans/sample_17


LC_ALL=C fgrep -w -f  $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/extra_elg_transc sample_17


# ZC416.8b.1	2071	1807.945	5.663982	103.000
# Y105E8A.7d.1	3028	2822.144	6.575818	186.663
# Y105E8A.7b.1	2946	2730.228	6.636242	182.243
# ZC416.8a.1	1920	1717.451	6.946512	120.000
# B0564.10b.1	1231	956.767	3.325168	32.000
# B0564.1a.2	1194	887.651	21.434201	191.372
# B0564.11.2	1284	909.639	16.722163	153.000
# B0564.1b.6	1207	902.121	41.396393	375.628

#ok, repair the "missing transcript" thing



cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/elegans/
rm 00_elg_trans_id
for i in *; do LC_ALL=C fgrep -w -v -f  $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/extra_elg_transc $i > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/02_grep_remove_extra_transc_elegans/$i; done &

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/02_grep_remove_extra_transc_elegans/

for i in *; do wc -l $i; done
	#10718 , woohoo!


#now, inopinata....

cd $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $3,$1}' 13_merge_og > 16_merge_og


mkdir $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/03_merge_inopinata

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/01_grep/inopinata

for i in *; do perl /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl -k -e "no_key" $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/16_merge_og $i  2> $i.merge_pl.error > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/03_merge_inopinata/$i; done &


for i in *; do grep "no_key" $i | wc  -l; done
	#0! woohoo!
	#also no errors, great

#


mkdir $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/04_awk_inopinata

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/03_merge_inopinata/

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6}' $i > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/04_awk_inopinata/$i; done

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/02_grep_remove_extra_transc_elegans


for i in *; do cp $i $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/04_awk_inopinata//$i; done

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/04_awk_inopinata


#cool

cd $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2}' 13_merge_og > 17_merge_og


mkdir $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/05_merge/


cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/04_awk_inopinata

for i in *; do perl /projects/phillipslab/gavincw/repeats_12-18-18/merge.pl -k -e "no_key" $wkdir/11_new_ortho_fork/02_get_1-1_inopinata_elegans_orthologs/17_merge_og $i  2> $i.merge_pl.error > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/05_merge/$i; done &


mkdir $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/06_awk

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/05_merge

for i in *; do awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2,$3,$4,$5,$6}' $i > $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/06_awk/$i; done

cd $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/06_awk/

for i in *; do grep "no_key" $i | wc -l; done 
	#0, woohoo! 


for i in *; do echo -e "Name\tLength\tEffectiveLength\tTPM\tNumReads" | cat - $i > $i.tmp && mv $i.tmp $i; done

#ok, now deseq


mkdir $wkdir/11_new_ortho_fork/04_DESeq2
mkdir $wkdir/11_new_ortho_fork/04_DESeq2/quants



while read F; do
  ln -s $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/06_awk/$F $wkdir/11_new_ortho_fork/04_DESeq2/quants/$F
done < $wkdir/sample_list.txt

ln -s $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/06_awk/sample_2 $wkdir/11_new_ortho_fork/04_DESeq2/quants/sample_2

ln -s $wkdir/11_new_ortho_fork/03_salmon_quant_1-1_orthologs/06_awk/sample_3 $wkdir/11_new_ortho_fork/04_DESeq2/quants/sample_3

grep -w -v sample_16 $wkdir/08_DESeq2/sample_table.tsv > $wkdir/11_new_ortho_fork/04_DESeq2/sample_table.tsv


awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1, $1}' $wkdir/11_new_ortho_fork/04_DESeq2/quants/sample_1 > $wkdir/11_new_ortho_fork/04_DESeq2/tx2gene.tsv


sed -i '1d' tx2gene.tsv
echo -e "TXNAME\tGENEID" | cat - tx2gene.tsv > tx2gene.tsv.tmp && mv tx2gene.tsv.tmp tx2gene.tsv

#should be ready for deseq2, see deseq2.R
