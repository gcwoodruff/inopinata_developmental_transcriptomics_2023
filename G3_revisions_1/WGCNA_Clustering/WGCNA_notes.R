##This is the code for WGCNA, with all C. inopinata genes (not just single-copy orthologs), for Woodruff et al. "Widespread changes in gene expression accompany body size evolution in nematodes"

#this workflow is based on a previous post: https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html

#load libraries
library(magrittr)
library(WGCNA)
library(ggplot2)


#set working directory

setwd("/Users/gavin/genome/transcriptome_revisions_12-2023/WGCNA/")

#get data in there...

	#transformed expression data; generated in DESeq2_all_genes_notes.R , https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/DESeq2_all_genes_notes.R
inop_tx_dat <- read.table("inopinata_all_genes_rlog_transformed_transcript_counts.tsv", sep="\t", header=TRUE,row.names=1)

	#inopinata gene ortholog categories; generated in join_orthogroups_genes_categories_notes.R , https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/join_orthogroups_genes_categories_notes.R

inop_ortho <- read.table("inopinata_GENES_ogs_categories.tsv", sep="\t", header=FALSE)
	
	#inopinata gene-domain info; generated in interproscan/workflow_notes.sh , https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/interproscan/workflow_notes.sh

inopinata_interpro <- read.table("inopinata_04.tsv", sep="\t", header=FALSE)

#transpose the data for the pickSoftThreshold() to pick an appropriate soft-thresholding power for network construction
tdat = t(inop_tx_dat)


#Analysis of scale free topology for picking an appropriate soft-thresholding power for network construction (see https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html) 

sft <- pickSoftThreshold(tdat,
  dataIsExpr = TRUE,
  corFnc = cor,
  networkType = "signed"
)

#visualize the sft output for choosing a soft-thresholding power parameter

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +  geom_point() +  geom_text(nudge_y = 0.1) +  geom_hline(yintercept = 0.80, col = "red") + ylim(c(min(sft_df$model_fit), 1.05)) + xlab("Soft Threshold (power)") + ylab("Scale Free Topology Model Fit, signed R^2") + theme_classic()

#by eyeballing it, let's do 12

#network construction (see https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html) 

bwnet <- blockwiseModules(tdat,
  maxBlockSize = 5000,
  TOMType = "signed", 
  power = 12, # soft threshold for network construction chosen via visualization above
  numericLabels = TRUE, 
  randomSeed = 1234, 
)
	#this may take a while

#write results so we don't have to construct the network again

readr::write_rds(bwnet,file =  "1-4-24_inopinata_wgcna_results.RDS")

#if we want to load the network again later, use the line below
#bwnet <- readRDS("1-4-24_inopinata_wgcna_results.RDS")

#get the modules
module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

#load data connected to the samples, inopinata_sample_table.tsv, can be found here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/inopinata_sample_table.tsv
metadata <- read.table("inopinata_sample_table.tsv",header=TRUE,sep='\t')

#get stage factor levels right
metadata$Stage <- factor(metadata$Stage,levels = c("L3", "L4", "adult"))

#does sample data align with network

all.equal(metadata$Sample_ID, rownames(module_eigengenes))
#[1] TRUE
#good
 
# Create the design matrix for linear models, we're interested in a relationship with developmental stage
des_mat <- model.matrix(~ metadata$Stage)

#Run linear model on each module. 
# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)


# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

#get modules with adjusted p-value < than 0.05
stats_df[stats_df$adj.P.Val < 0.05,]
	#these are the modules we look at further

#how many are there
nrow(stats_df[stats_df$adj.P.Val < 0.05,])
#[1] 12

#let's write these results to a file , this is supplemental table 16
write.table(stats_df[order(-stats_df$metadata.Stageadult,stats_df$adj.P.Val),], "WGCNA_cluster_linear_model_coefficients.tsv",quote=FALSE,sep="\t")

#get the significant modules

sig_modules <- stats_df[stats_df$adj.P.Val < 0.05,]$module


#ok, merge modules to genes

inop_tx_dat$gen_id <- row.names(inop_tx_dat)

gene_module_key <- tibble::enframe(bwnet$colors, name = "gen_id", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

#write the constituents of the modules to a file
write.table(gene_module_key, "WGCNA_modules_genes.tsv",quote=FALSE,sep="\t")
  #this is supplemental table 13

#connect transcript abundances to the modules
inop_tx_dat_modules <- merge(inop_tx_dat,gene_module_key)


#ok, melt the data for plotting with ggplot
library(reshape2)
meltdat <- melt(inop_tx_dat_modules,id.var=c("gen_id","module"))

#get column names right
names(meltdat)[names(meltdat) == 'variable'] <- 'Sample_ID'
names(meltdat)[names(meltdat) == 'value'] <- 'expression'

#combine this with the sample data
meltdat_merge <- merge(meltdat,metadata)


#like with supplemental figure 3A, I want to make line graphs showing the average gene expression at each stage. To do this, we need to get averages among samples at the same stage. This is what I'm doing here

#make a column with all the info we need-- stage, gene id, and module id ; we'll aggregate by this column
meltdat_merge$stage.gene.module <- paste(meltdat_merge$Stage,meltdat_merge$gen_id,meltdat_merge$module)

#get the average per stage per gene (and we also have the module id)
meltdat_merge_agg <- aggregate(expression ~ stage.gene.module, FUN=mean,data=meltdat_merge)

#now, split up the columns so we can plot it!
library(tidyr)
meltdat_merge_agg<- separate(meltdat_merge_agg, stage.gene.module, into = c("stage.gene","module"), sep = " (?=[^ ]+$)")

meltdat_merge_agg<- separate(meltdat_merge_agg, stage.gene, into = c("stage","gene"), sep = " (?=[^ ]+$)")

#let's plot some things

#load some libraries i like for this
library(cowplot)
library(lemon)

#get stage factor levels right
meltdat_merge_agg$stage <- factor(meltdat_merge_agg$stage,levels = c("L3", "L4", "adult"))


#this plots ALL the modules. I did not include this in the paper, but it's here.
ggplot(meltdat_merge_agg, aes(x=stage, y=expression,group=gene)) + geom_line(alpha=0.15,size=0.15) + facet_rep_wrap(~module) + scale_color_brewer(palette="Set1") + theme_cowplot()



#let's just look at modules with a significant association with developmental change


meltdat_merge_agg_sig <- droplevels(  meltdat_merge_agg[meltdat_merge_agg$module %in% sig_modules,])

#get module factors right. that is, in an order that makes sense
meltdat_merge_agg_sig$module <- factor(meltdat_merge_agg_sig$module,levels = c("ME1","ME2","ME3","ME4","ME7","ME10","ME11","ME14","ME15","ME26","ME37","ME38"))

#this is supplemental figure 4
ggplot(meltdat_merge_agg_sig, aes(x=stage, y=expression,group=gene)) + geom_line(alpha=0.15,size=0.15) + facet_rep_wrap(~module) + scale_color_brewer(palette="Set1") + theme_cowplot()

#okay, now let's look at ortholog categories ; this is very similar for what was done for genes with positive and negative trajectories over developmental time. making supplemental figure 5 and supplemental table 15

#let's give the gene column an the appropriate name for the OG and TX df's so we can merge them
inop_ortho$gene <- inop_ortho$V1

inop_tx_dat_modules$gene <- inop_tx_dat_modules$gen_id

#merge the df's
ortho_merge <- merge(inop_tx_dat_modules,inop_ortho, by.x="gene", by.y="gene")

#get total number of ortholog category counts across whole genome
total_ortho_cat <- as.data.frame(table(inop_ortho$V3))
total_ortho_cat$module <- "All genes"

#split the df by the significant modules
ME02 <- ortho_merge[ortho_merge$module == "ME2",]
ME11 <- ortho_merge[ortho_merge$module == "ME11",]
ME01<- ortho_merge[ortho_merge$module=="ME1",]
ME04 <- ortho_merge[ortho_merge$module == "ME4",]
ME10 <- ortho_merge[ortho_merge$module == "ME10",]
ME03 <- ortho_merge[ortho_merge$module == "ME3",]
ME14 <- ortho_merge[ortho_merge$module == "ME14",]
ME07 <- ortho_merge[ortho_merge$module == "ME7",]
ME38 <- ortho_merge[ortho_merge$module == "ME38",]
ME15 <- ortho_merge[ortho_merge$module == "ME15",]
ME37 <- ortho_merge[ortho_merge$module == "ME37",]
ME26 <- ortho_merge[ortho_merge$module == "ME26",]

#count the number of genes in ortholog category in each module
ME02_ortho_cat <- as.data.frame(table(ME02$V3))
ME11_ortho_cat <- as.data.frame(table(ME11$V3))
ME01_ortho_cat <- as.data.frame(table(ME01$V3))
ME04_ortho_cat <- as.data.frame(table(ME04$V3))
ME10_ortho_cat <- as.data.frame(table(ME10$V3))
ME03_ortho_cat <- as.data.frame(table(ME03$V3))
ME14_ortho_cat <- as.data.frame(table(ME14$V3))
ME07_ortho_cat <- as.data.frame(table(ME07$V3))
ME38_ortho_cat <- as.data.frame(table(ME38$V3))
ME15_ortho_cat <- as.data.frame(table(ME15$V3))
ME37_ortho_cat <- as.data.frame(table(ME37$V3))
ME26_ortho_cat <- as.data.frame(table(ME26$V3))

#give the df the module id
ME02_ortho_cat$module <-  "ME2"
ME11_ortho_cat$module <- "ME11"
ME01_ortho_cat$module <- "ME1"
ME04_ortho_cat$module <-  "ME4"
ME10_ortho_cat$module <- "ME10"
ME03_ortho_cat$module <-  "ME3"
ME14_ortho_cat$module <- "ME14"
ME07_ortho_cat$module <-  "ME7"
ME38_ortho_cat$module <- "ME38"
ME15_ortho_cat$module <- "ME15"
ME37_ortho_cat$module <- "ME37"
ME26_ortho_cat$module <- "ME26"

#combine the ortholog category count df's

ortho_cat_df <- rbind(total_ortho_cat,ME02_ortho_cat,ME11_ortho_cat,ME01_ortho_cat,ME04_ortho_cat,ME10_ortho_cat,ME03_ortho_cat,ME14_ortho_cat,ME07_ortho_cat,ME38_ortho_cat,ME15_ortho_cat,ME37_ortho_cat,ME26_ortho_cat)

#get module factor right
ortho_cat_df$module <- as.factor(ortho_cat_df$module)

ortho_cat_df$module <- factor(ortho_cat_df$module,levels = c("All genes","ME1","ME2","ME3","ME4","ME7","ME10","ME11","ME14","ME15","ME26","ME37","ME38"))

#count the total number of genes in each relevant module (including whole genome)
tot_genes <- aggregate(Freq ~ module, FUN=sum, data=ortho_cat_df)


#plot this!


ggplot(ortho_cat_df, aes(x = module, y = Freq, fill = Var1)) + geom_col(position = "fill") + scale_fill_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="bottom",legend.title=element_blank(),legend.justification = "center") + annotate("text",x=1, y=0.05, label=tot_genes$Freq[1],size=5) + annotate("text",x=2, y=0.05, label=tot_genes$Freq[2],size=5) + annotate("text",x=3, y=0.05, label=tot_genes$Freq[3],size=5) + annotate("text",x=4, y=0.05, label=tot_genes$Freq[4],size=5) + annotate("text",x=5, y=0.05, label=tot_genes$Freq[5],size=5) + annotate("text",x=6, y=0.05, label=tot_genes$Freq[6],size=5) + annotate("text",x=7, y=0.05, label=tot_genes$Freq[7],size=5) + annotate("text",x=8, y=0.05, label=tot_genes$Freq[8],size=5)+ annotate("text",x=9, y=0.05, label=tot_genes$Freq[9],size=5)+ annotate("text",x=10, y=0.05, label=tot_genes$Freq[10],size=5)+ annotate("text",x=11, y=0.05, label=tot_genes$Freq[11],size=5)+ annotate("text",x=12, y=0.05, label=tot_genes$Freq[12],size=5)+ annotate("text",x=13, y=0.05, label=tot_genes$Freq[13],size=5) + xlab("Module") + ylab("Proportion") 

#this is supplemental figure 5!

#next, chi square tests (chi-square tests)

#load a library we need
library(dplyr)
#cast the data by module to get it in a form we want for this purpose
orthocast <- dcast(ortho_cat_df,Var1~module,value.var="Freq")
#replace na with 0
orthocast[is.na(orthocast)] <- 0


#make an empty variable that we'll add stuff (chi square output for all modules) to for the loop
testfun<-NULL

#a loop and a loop within that loop

	#for all modules
for(i in levels(ortho_cat_df$module)[-1]){
	#this is the module of interest
  the_i <- i
  #make an empty variable that we'll add stuff (ie chi square stats) the inner loop
  i_df <-NULL
  #make the df we need for doing chi square tests...
	#the ortholog category
	#the number of genes in the genome in those categories
	#the number of genes in the genome NOT in those categories
	#the number of genes in the module of interest in those categories
	#the number of genes in the module of interest NOT in those categories
  i_df <- data.frame(orthotype=orthocast$Var1, All_genes=orthocast[,2],All_genes_no= tot_genes[tot_genes$module == "All genes",]$Freq-orthocast[,2],positive=orthocast %>% select(all_of(the_i)),positive_no=tot_genes[tot_genes$module == the_i,]$Freq-( orthocast %>% select(all_of(the_i))))

  #now another loop! here, for all ortholog categories in that module...
  for(j in 1:nrow(i_df)){
  	#the ortholog category row
    orthoclass <- i_df$orthotype[j]
    #make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that are and are not that specific ortholog category in that specific module compared with all genes in the genome
    the_matrix <- matrix(c(i_df[j,2],i_df[j,3],i_df[j,4],i_df[j,5]),ncol=2)
    #do a chi square test on that matrix
    matrixtest <- chisq.test(the_matrix)
    #make a row for our df including:
		#the orthogroup category
		#the number of genes in the genome that are in that category
		#the number of genes in the genome that are not in that category
		#the number of module-specific genes in that category
		#the number of module-specific trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of module-specific genes in that category
		#the chi square test statistic
		#the chisquare p value
		#the module of interest
    row_to_add <- data.frame(orthotype=orthoclass,all_genes_with=i_df[j,2],all_genes_without=i_df[j,3], test_genes_with=i_df[j,4],test_genes_without=i_df[j,5], fra_all_genes_with=i_df[j,2]/sum(i_df[,2]), fra_positive_genes_with=i_df[j,4]/sum(i_df[,4]), chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value,module=the_i)
    testfun <- rbind(testfun,row_to_add)
  }

}

#adjust p values
testfun$p.adjust <- p.adjust(testfun$chisq_p,method="BH")
#write this table to a file
write.table(testfun,"WGCNA_modules_ortholog_type_chi_square.tsv",quote=FALSE,row.names = FALSE,sep="\t")
	#this is supplemental table 15!

#now, domains

#get domain data column names right
names(inopinata_interpro)[names(inopinata_interpro) == 'V1'] <- 'gene'
names(inopinata_interpro)[names(inopinata_interpro) == 'V2'] <- 'interpro_id'
names(inopinata_interpro)[names(inopinata_interpro) == 'V3'] <- 'interpro_description'

#merge domain data with expression data

big_df_merge <- merge(inopinata_interpro,ortho_merge)

#count the domains

inopinata_domain_counts <- as.data.frame(table(inopinata_interpro$interpro_id))

#get a domain description df
interpro_domains_descriptions <- data.frame(interpro_id=inopinata_interpro$interpro_id,interpro_description=inopinata_interpro$interpro_description)
nrow(interpro_domains_descriptions)
#remove duplicate entries
interpro_domains_descriptions <- interpro_domains_descriptions[!duplicated(interpro_domains_descriptions), ]
nrow(interpro_domains_descriptions)


#ok, join descriptions and get some frequencies, differences
#get column names right for the domain count table
names(inopinata_domain_counts)[names(inopinata_domain_counts) == 'Var1'] <- 'interpro_id'
names(inopinata_domain_counts)[names(inopinata_domain_counts) == 'Freq'] <- 'count'

#merge the counts with the interpro domain descriptions

inopinata_domain_counts <- merge(inopinata_domain_counts,interpro_domains_descriptions)

#get total number of gene-domain pairs
inopinata_domain_counts$total <- sum(inopinata_domain_counts$count)
#get gene-domain pair counts that are NOT of a given domain
inopinata_domain_counts$num_not_domain <- inopinata_domain_counts$total-inopinata_domain_counts$count
#get the fraction of gene-domain counts that are in a given domain
inopinata_domain_counts$fra_domain <- inopinata_domain_counts$count/inopinata_domain_counts$total



#okay, do what we did before (with OG categories), but this time with domains-- split by module, looking only at modules with significant relationships with developmental time
	#here, the whole genome
inopinata_domain_counts$module <- "All genes"
#here, all sig modules
ME02 <- big_df_merge[big_df_merge$module == "ME2",]
ME11 <- big_df_merge[big_df_merge$module == "ME11",]
ME01<- big_df_merge[big_df_merge$module=="ME1",]
ME04 <- big_df_merge[big_df_merge$module == "ME4",]
ME10 <- big_df_merge[big_df_merge$module == "ME10",]
ME03 <- big_df_merge[big_df_merge$module == "ME3",]
ME14 <- big_df_merge[big_df_merge$module == "ME14",]
ME07 <- big_df_merge[big_df_merge$module == "ME7",]
ME38 <- big_df_merge[big_df_merge$module == "ME38",]
ME15 <- big_df_merge[big_df_merge$module == "ME15",]
ME37 <- big_df_merge[big_df_merge$module == "ME37",]
ME26 <- big_df_merge[big_df_merge$module == "ME26",]

#count the domains among the genes in the significant modules


ME02_ortho_cat <- as.data.frame(table(ME02$interpro_id))
ME11_ortho_cat <- as.data.frame(table(ME11$interpro_id))
ME01_ortho_cat <- as.data.frame(table(ME01$interpro_id))
ME04_ortho_cat <- as.data.frame(table(ME04$interpro_id))
ME10_ortho_cat <- as.data.frame(table(ME10$interpro_id))
ME03_ortho_cat <- as.data.frame(table(ME03$interpro_id))
ME14_ortho_cat <- as.data.frame(table(ME14$interpro_id))
ME07_ortho_cat <- as.data.frame(table(ME07$interpro_id))
ME38_ortho_cat <- as.data.frame(table(ME38$interpro_id))
ME15_ortho_cat <- as.data.frame(table(ME15$interpro_id))
ME37_ortho_cat <- as.data.frame(table(ME37$interpro_id))
ME26_ortho_cat <- as.data.frame(table(ME26$interpro_id))

#give the dfs module ids


ME02_ortho_cat$module <-  "ME2"
ME11_ortho_cat$module <- "ME11"
ME01_ortho_cat$module <- "ME1"
ME04_ortho_cat$module <-  "ME4"
ME10_ortho_cat$module <- "ME10"
ME03_ortho_cat$module <-  "ME3"
ME14_ortho_cat$module <- "ME14"
ME07_ortho_cat$module <-  "ME7"
ME38_ortho_cat$module <- "ME38"
ME15_ortho_cat$module <- "ME15"
ME37_ortho_cat$module <- "ME37"
ME26_ortho_cat$module <- "ME26"

#combine the dfs

domain_cat_df <- rbind(ME02_ortho_cat,ME11_ortho_cat,ME01_ortho_cat,ME04_ortho_cat,ME10_ortho_cat,ME03_ortho_cat,ME14_ortho_cat,ME07_ortho_cat,ME38_ortho_cat,ME15_ortho_cat,ME37_ortho_cat,ME26_ortho_cat)


#chi square stats, domains

inopinata_domain_counts_ii <- as.data.frame(table(inopinata_interpro$interpro_id))
inopinata_domain_counts_ii$module <- "All_genes"
#join module domain info with whole genome domain info
domain_cat_df <- rbind(domain_cat_df,inopinata_domain_counts_ii)
#count the total number of gene-domain pairs in each module (including whole genome)
tot_genes <- aggregate(Freq ~ module, FUN=sum, data=domain_cat_df)

#cast the data to get it in a form we need for chi square testing
orthocast <- dcast(domain_cat_df,Var1~module,value.var="Freq")
#replace NA with 0
orthocast[is.na(orthocast)] <- 0
#get factors right
domain_cat_df$module <- as.factor(domain_cat_df$module)


#make an empty variable that we'll add stuff (chi square output for all modules) to for the loop
testfun<-NULL
#a loop and a loop within that loop

	#for all modules

for(i in levels(domain_cat_df$module)[-1]){
	#this is the module of interest
  the_i <- i
  #make an empty variable that we'll add stuff (ie chi square stats) the inner loop
  i_df <-NULL
  #make the df we need for doing chi square tests...
	#the domain
	#the number of genes in the genome with that domain
	#the number of genes in the genome that do NOT have that domain
	#the number of genes in the module of interest with that domain
	#the number of genes in the module of interest that do NOT have that domain
  i_df <- data.frame(orthotype=orthocast$Var1, All_genes=orthocast[,2],All_genes_no= tot_genes[tot_genes$module == "All_genes",]$Freq-orthocast[,2],positive=orthocast %>% select(all_of(the_i)),positive_no=tot_genes[tot_genes$module == the_i,]$Freq-(orthocast %>% select(all_of(the_i))))
  #now another loop! here, for all domains in that specific module...
  for(j in 1:nrow(i_df)){
  	#the domain (row)
    orthoclass <- i_df$orthotype[j]
    #make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that do and do not have that specific domain in that specific module compared with all genes in the genome
    the_matrix <- matrix(c(i_df[j,2],i_df[j,3],i_df[j,4],i_df[j,5]),ncol=2)
    #do a chi square test on that matrix
    matrixtest <- chisq.test(the_matrix)
    #make a row for our df including:
		#the domain
		#the number of genes in the genome with that domain
		#the number of genes in the genome that do not have that domain
		#the number of module-specific genes with that domain
		#the number of module-specific genes that do not have that domain
		#the fraction of genes in the genome with that domain (among all gene-domain pairs!)
		#the fraction of module-specific genes with that domain (among all gene-domain pairs that module has!)
		#the chi square test statistic
		#the chisquare p value
		#the module of interest
    row_to_add <- data.frame(orthotype=orthoclass,all_genes_with=i_df[j,2],all_genes_without=i_df[j,3], test_genes_with=i_df[j,4],test_genes_without=i_df[j,5], fra_all_genes_with=i_df[j,2]/sum(i_df[,2]), fra_positive_genes_with=i_df[j,4]/sum(i_df[,4]), chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value,module=the_i)
    testfun <- rbind(testfun,row_to_add)
  }

}

	#this may take a while

#adjust p value
testfun$p.adjust <- p.adjust(testfun$chisq_p,method="BH")
#as i reused a loop, change orthotype to interpro_id
names(testfun)[names(testfun) == 'orthotype'] <- 'interpro_id'

#merge with useful descriptions
testfun_merge <- merge(testfun,interpro_domains_descriptions)
#write to a table
#not in the paper, but is on the github
write.table(testfun_merge,"WGCNA_modules_domains_chi_square_with_descriptions.tsv",quote=FALSE,row.names = FALSE,sep="\t")

#get domains with significant enrichment/depletion across modules

sig_ids <- testfun_merge[testfun_merge$p.adjust < 0.05,]
write.table(sig_ids,"WGCNA_sig_interpro_ids_with_descriptions.tsv",quote=FALSE,row.names = FALSE,sep="\t")
	#this is supplemental table 14










