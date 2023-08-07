
#interaction list
my_RNAseq_genes <- scan("inopinata_interaction_list.txt",what = character(),quote="")

liang_genes <- scan("liang_list_sma-9_and_dbl-1_repair_2021.txt",what = character(),quote="")


roberts_genes <- scan("roberts_2010_bmc_dev_bio_list.txt",what = character(),quote="")


lakdawala_genes <- scan("Lakdawala_2018_MBoC_list.txt",what = character(),quote="")


interaction_RNAseq_genes<- unique(my_RNAseq_genes)
liang_genes<- unique(liang_genes)
roberts_genes<-unique(roberts_genes)
lakdawala_genes <- unique(lakdawala_genes)

num_interaction_genes <- length(interaction_RNAseq_genes)
num_liang_genes <- length(liang_genes)
num_roberts_genes<- length(roberts_genes)
num_lakdawala_genes <- length(lakdawala_genes)
num_genes_genome <- 20040

liang_num_intersect <- length(intersect(my_RNAseq_genes,liang_genes))
roberts_num_intersect <- length(intersect(my_RNAseq_genes,roberts_genes))
lakdawala_num_intersect <- length(intersect(my_RNAseq_genes, lakdawala_genes))

#test for liang list

phyper(liang_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_liang_genes,lower.tail = FALSE, log.p = FALSE)

# [1] 0.0910099
#p=0.09

#test for roberts list

phyper(roberts_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_roberts_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.9999788


#test for lakdawala list

phyper(lakdawala_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_lakdawala_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.1376305

#just L3 DFE elegans and inopinata

my_RNAseq_genes <- scan("L3_list.txt",what = character(),quote="")


interaction_RNAseq_genes<- unique(my_RNAseq_genes)

num_interaction_genes <- length(interaction_RNAseq_genes)
liang_num_intersect <- length(intersect(my_RNAseq_genes,liang_genes))
roberts_num_intersect <- length(intersect(my_RNAseq_genes,roberts_genes))
lakdawala_num_intersect <- length(intersect(my_RNAseq_genes, lakdawala_genes))
#test for liang list
phyper(liang_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_liang_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.5429832

#test for roberts list

phyper(roberts_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_roberts_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.9999225


#test for lakdawala list

phyper(lakdawala_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_lakdawala_genes,lower.tail = FALSE, log.p = FALSE)
#[1] 0.4268777


#just L4 DFE elegans and inopinata

my_RNAseq_genes <- scan("L4_list.txt",what = character(),quote="")


interaction_RNAseq_genes<- unique(my_RNAseq_genes)

num_interaction_genes <- length(interaction_RNAseq_genes)
liang_num_intersect <- length(intersect(my_RNAseq_genes,liang_genes))
roberts_num_intersect <- length(intersect(my_RNAseq_genes,roberts_genes))
lakdawala_num_intersect <- length(intersect(my_RNAseq_genes, lakdawala_genes))
#test for liang list
phyper(liang_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_liang_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.1088713

#test for roberts list

phyper(roberts_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_roberts_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.9999905


#test for lakdawala list

phyper(lakdawala_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_lakdawala_genes,lower.tail = FALSE, log.p = FALSE)
#[1] 0.01503604

#just Adult DFE elegans and inopinata

my_RNAseq_genes <- scan("adult_list.txt",what = character(),quote="")


interaction_RNAseq_genes<- unique(my_RNAseq_genes)

num_interaction_genes <- length(interaction_RNAseq_genes)
liang_num_intersect <- length(intersect(my_RNAseq_genes,liang_genes))
roberts_num_intersect <- length(intersect(my_RNAseq_genes,roberts_genes))
lakdawala_num_intersect <- length(intersect(my_RNAseq_genes, lakdawala_genes))
#test for liang list
phyper(liang_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_liang_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.04815952

#test for roberts list

phyper(roberts_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_roberts_genes,lower.tail = FALSE, log.p = FALSE)

#[1] 0.9999931


#test for lakdawala list

phyper(lakdawala_num_intersect,num_interaction_genes,num_genes_genome-num_interaction_genes,num_lakdawala_genes,lower.tail = FALSE, log.p = FALSE)
#[1] 0.4806446

#ok, adjust some p-values

p_vals <- p.adjust(sort(c(0.0910099,0.9999788,0.1376305,0.5429832,0.9999225,0.4268777,0.1088713,0.9999905,0.01503604,0.04815952,0.9999931,0.4806446)))

# [1] 0.1804325 0.5297547 0.9100990 0.9798417 1.0000000 1.0000000 1.0000000
# [8] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000

p_vals <- p.adjust(sort(c(0.0910099,0.9999788,0.1376305,0.5429832,0.9999225,0.4268777,0.1088713,0.9999905,0.01503604,0.04815952,0.9999931,0.4806446)), method="BH")
#> p_vals
# [1] 0.1804325 0.2889571 0.3266139 0.3266139 0.3303132 0.8144748 0.8144748
# [8] 0.8144748 0.9999931 0.9999931 0.9999931 0.9999931

#q=size of overlap-1;
#m=number of upregulated genes in experiment #1;
#n=(total number of genes on platform-m);
#k=number of upregulated genes in experiment #2. 
#


#q=size of overlap-1;
	#9-1 = 8
#m=number of upregulated genes in experiment #1;
	#201
#n=(total number of genes on platform-m);
	# 20040 - 201 = 19839
#k=number of upregulated genes in experiment #2. 
	#2438
#

#phyper(8, 201, 19839, 2438, lower.tail = FALSE, log.p = FALSE)
#
#	# 0.9999517 ; no significant overlap
#
##what if all overlap??
#phyper(200, 201, 19839, 2438, lower.tail = FALSE, log.p = FALSE)
#	# [1] 7.352967e-188 , cool
#
#
##doing this again ; here, liang list and my 10% list
#
##q=size of overlap-1;
#	#9-1 = 8
##m=number of upregulated genes in experiment #1;
#	#440
##n=(total number of genes on platform-m);
#	# 10717 - 440 = 10277
##k=number of upregulated genes in experiment #2. 
#	#272
##
#
#
#phyper(8, 440, 10277, 272, lower.tail = FALSE, log.p = FALSE)
#
#	#[1] 0.7915795 , damn still no significant overlap
#
#
#
##doing this again ; here, roberts list and my 10% list
#
##q=size of overlap-1;
#	#22-1 = 21
##m=number of upregulated genes in experiment #1;
#	#440
##n=(total number of genes on platform-m);
#	# 10717 - 440 = 10277
##k=number of upregulated genes in experiment #2. 
#	#2436
##
#
#
#phyper(21, 440, 10277, 2436, lower.tail = FALSE, log.p = FALSE)
#
#	#[1] 1 , still no significant overlap
