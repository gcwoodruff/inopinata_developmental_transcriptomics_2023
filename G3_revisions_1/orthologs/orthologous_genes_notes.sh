#This is how orthogroups were categorized for Woodruff et al. "Widespread changes in gene expression accompany body size evolution in nematodes"

#This is the first part, the next part is "join_orthogroups_genes_categories_notes.R"

#the previous orthofinder run for this is here: /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/inopinata_RNA-seq_phillipslab/inopinata_RNAseq_1-20-20/11_new_ortho_fork/01_get_largest_isoform/OrthoFinder/Results_Feb20/WorkingDirectory/OrthoFinder/Results_Feb20

#gene count file is at Orthogroups.GeneCount.tsv

#columns for the orthogroup files
#1 Orthogroup
#2 briggsae
#3 elegans
#4 inopinata
#5 nigoni
#6 remanei
#7 Total

#set working directory

wkdir="/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/transcriptomics_revisions_1-2024"

#okay, let's organize some types of orthogroups

#the single copy-orthogroups between elegans and inopinata...

mkdir $wkdir/01_orthogroups

mkdir $wkdir/01_orthogroups/00_awk_single-copy/

#get the elegans-inopinata single-copy orthogroups

#Orthogroups.GeneCount.tsv is here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/Orthogroups.GeneCount.tsv

awk 'BEGIN {FS="\t"} $3 == "1" && $4 == "1" {print} ' Orthogroups.GeneCount.tsv > $wkdir/01_orthogroups/00_awk_single-copy/single-copy_orthologs_eleg_inop.tsv

#how many

cd $wkdir/01_orthogroups/00_awk_single-copy/

wc -l single-copy_orthologs_eleg_inop.tsv
#10719 single-copy_orthologs_eleg_inop.tsv
#sweet, that worked. that's the number from the initial submission. okay, let's get some more genes.

mkdir $wkdir/01_orthogroups/01_awk_multi-copy/

#orthogroups where elegans and inopinata have more than one copy
awk 'BEGIN {FS="\t"} $3 > 1 && $4 > 1 {print} ' Orthogroups.GeneCount.tsv > $wkdir/01_orthogroups/01_awk_multi-copy/multi-copy_orthologs_both_eleg_inop.tsv

#orthogroups where elegans has one copy but inopinata has more than one copy
awk 'BEGIN {FS="\t"} $3 == 1 && $4 > 1 {print} ' Orthogroups.GeneCount.tsv > $wkdir/01_orthogroups/01_awk_multi-copy/multi-copy_inop_single-copy_eleg.tsv

#orthogroups where elegans has more than one copy but inopinata has one copy
awk 'BEGIN {FS="\t"} $3 > 1 && $4 == 1 {print} ' Orthogroups.GeneCount.tsv > $wkdir/01_orthogroups/01_awk_multi-copy/multi-copy_eleg_single-copy_inop.tsv

#using venny to show these are completely non-overlapping sets. they are, great!
#awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' multi-copy_orthologs_both_eleg_inop.tsv
#awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' multi-copy_inop_single-copy_eleg.tsv
#awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' multi-copy_eleg_single-copy_inop.tsv

#get genes found in other caenorhabditis but not elegans or inopinata
mkdir $wkdir/01_orthogroups/02_awk_zero_copy_number/

#get orthogroups where elegans has zero copies, but inopinata has more than zero copies
awk 'BEGIN {FS="\t"} $3 == 0 && $4 > 0 {print} ' Orthogroups.GeneCount.tsv > $wkdir/01_orthogroups/02_awk_zero_copy_number/elegans_zero_copy_inopinata_present.tsv

cd $wkdir/01_orthogroups/02_awk_zero_copy_number/
#get those orthogroups that are multicopy but present ONLY in inopinata
awk 'BEGIN {FS="\t"} $2 ==0 && $3 == 0 && $5 == 0 && $6 == 0 {print} ' elegans_zero_copy_inopinata_present.tsv > inopinata_specific_multicopy.tsv

#get list of inopinata-specific multicopy ogs
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' inopinata_specific_multicopy.tsv > inopinata_specific_multicopy_OGS.txt
#use grep to extract orthogroups where elegans has zero copies, inopinata has more than zero copies, but another caeno species also has at least one copy
grep -w -v -f inopinata_specific_multicopy_OGS.txt elegans_zero_copy_inopinata_present.tsv > multicaeno_elegans_zero_copy_inopinata_present.tsv

#do the numbers add up

wc -l elegans_zero_copy_inopinata_present.tsv
#910
wc -l inopinata_specific_multicopy.tsv
#364
wc -l multicaeno_elegans_zero_copy_inopinata_present.tsv
#546
#> 364+546
#[1] 910
#they add up i think this worked. now, let's do stuff not present in inopinata

awk 'BEGIN {FS="\t"} $3 > 0 && $4 == 0 {print} ' Orthogroups.GeneCount.tsv > $wkdir/01_orthogroups/02_awk_zero_copy_number/inopinata_zero_copy_elegans_present.tsv

#get those orthogroups that are multicopy but present ONLY in elegans
awk 'BEGIN {FS="\t"} $2 ==0 && $4 == 0 && $5 == 0 && $6 == 0 {print} ' inopinata_zero_copy_elegans_present.tsv > elegans_specific_multicopy.tsv

#get list of elegans-specific multicopy ogs
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' elegans_specific_multicopy.tsv > elegans_specific_multicopy_OGS.txt

#use grep to extract orthogroups where inopinata has zero copies, elegans has more than zero copies, but another caeno species also has at least one copy
grep -w -v -f elegans_specific_multicopy_OGS.txt inopinata_zero_copy_elegans_present.tsv > multicaeno_inopinata_zero_copy_elegans_present.tsv

wc -l inopinata_zero_copy_elegans_present.tsv
#1848
wc -l elegans_specific_multicopy.tsv
#307
wc -l multicaeno_inopinata_zero_copy_elegans_present.tsv
#1541
#> 307+1541
#[1] 1848
#okay it adds up

#now, completely orphan genes

mkdir $wkdir/01_orthogroups/03_orphan_genes/

#Orthogroups_UnassignedGenes.tsv is here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/Orthogroups_UnassignedGenes.tsv

#elegans orphan genes
awk 'BEGIN {FS="\t"} $3 != "" {print}' Orthogroups_UnassignedGenes.tsv > $wkdir/01_orthogroups/03_orphan_genes/elegans_orhpan_genes.tsv

#inopinata orphan genes
awk 'BEGIN {FS="\t"} $4 != "" {print}' Orthogroups_UnassignedGenes.tsv > $wkdir/01_orthogroups/03_orphan_genes/inopinata_orhpan_genes.tsv

#okay, categories of orthogroup
# "single-copy ortholog"
# "multi-copy, present in both species"
# "multi-caeno inopinata-absent"
# "multi-caeno elegans-absent"
# "inopinata-specific multi-copy"
# "elegans-specific multi-copy"
# "elegans orphan"
# "inopinata orphan"

mkdir $wkdir/01_orthogroups/04_add_copy-number_categories/
#add the category label "single-copy ortholog"
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"single-copy ortholog"}' $wkdir/01_orthogroups/00_awk_single-copy/single-copy_orthologs_eleg_inop.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/single_copy_orthogroups.tsv

cd $wkdir/01_orthogroups/01_awk_multi-copy/
#combine the multi-copy og files
	#remove header
sed -i '1d' multi-copy_orthologs_both_eleg_inop.tsv
#combine the multi-copy og files
cat * > multi-copy_orthologs.tsv
#add labels to orthogroup ID's
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"multi-copy orthogroup present in both species"}' $wkdir/01_orthogroups/01_awk_multi-copy/multi-copy_orthologs.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/multi_copy_orthogroups_present_in_both_species.tsv

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"multi-caeno inopinata-absent"}' $wkdir/01_orthogroups/02_awk_zero_copy_number/multicaeno_inopinata_zero_copy_elegans_present.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/multi-caeno_inopinata-absent_orthogroups.tsv

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"multi-caeno elegans-absent"}' $wkdir/01_orthogroups/02_awk_zero_copy_number/multicaeno_elegans_zero_copy_inopinata_present.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/multi-caeno_elegans-absent_orthogroups.tsv

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"inopinata-specific multi-copy"}' $wkdir/01_orthogroups/02_awk_zero_copy_number/inopinata_specific_multicopy.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/inopinata-specific_multi-copy_orthogroups.tsv

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"elegans-specific multi-copy"}' $wkdir/01_orthogroups/02_awk_zero_copy_number/elegans_specific_multicopy.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/elegans-specific_multi-copy_orthogroups.tsv

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"elegans orphan"}' $wkdir/01_orthogroups/03_orphan_genes/elegans_orhpan_genes.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/elegans_orphan_orthogroups.tsv

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,"inopinata orphan"}' $wkdir/01_orthogroups/03_orphan_genes/inopinata_orhpan_genes.tsv > $wkdir/01_orthogroups/04_add_copy-number_categories/inopinata_orphan_orthogroups.tsv

cd  $wkdir/01_orthogroups/04_add_copy-number_categories/
#remove headers
sed -i '1d' elegans_orphan_orthogroups.tsv
sed -i '1d' inopinata_orphan_orthogroups.tsv
#combine the files
cat * > orthogroup_categories.tsv
#test-- are there any duplicate orthogroups?
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1}' orthogroup_categories.tsv > orthogroups_test.txt
sort orthogroups_test.txt | uniq > orthogroups_test_sort_uniq.txt

wc -l orthogroups_test.txt
#17338 orthogroups_test.txt

wc -l orthogroups_test_sort_uniq.txt
#17338 orthogroups_test_sort_uniq.txt

#nice-- there are no duplicate orthogroups


#okay, getting everything in order here. a list of genes, orthogroups, and orthogroup categories.

cd  $wkdir/01_orthogroups/04_add_copy-number_categories/

#Orthogroups.tsv is here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/Orthogroups.tsv

#get the genes associated with my orthogroups
grep -w -f orthogroups_test_sort_uniq.txt Orthogroups.tsv > orthogroups_genes.tsv
#orphan genes
cd $wkdir/01_orthogroups/03_orphan_genes
#remove headers
sed -i '1d' elegans_orhpan_genes.tsv
sed -i '1d' inopinata_orhpan_genes.tsv
#combine orthogroups with more than one gene with orthogroups defined by orphan genes
cd $wkdir/01_orthogroups/04_add_copy-number_categories/

cat orthogroups_genes.tsv $wkdir/01_orthogroups/03_orphan_genes/elegans_orhpan_genes.tsv $wkdir/01_orthogroups/03_orphan_genes/inopinata_orhpan_genes.tsv > orthogroups_genes_w_orphans.tsv
#okay, has correct number of lines
#sort and get unique orthogroups for bash
sort orthogroups_genes_w_orphans.tsv | uniq > orthogroups_genes_w_orphans_sort_uniq.tsv
sort orthogroup_categories.tsv | uniq > orthogroup_categories_sort_uniq.tsv

#get all these files to local for R

###okay, connect these things-- finished this in R.










