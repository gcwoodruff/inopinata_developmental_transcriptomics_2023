#This is the code for addressing _all_ C. inopinata genes (not just single-copy orthologs) for Woodruff et al. "Widespread changes in gene expression accompany body size evolution in nematodes"

#here, DESeq2, linear models, plots for data that look at all genes by transcriptional dynamics category (positive, negative, no change across development).
	#additionally, chi-square to get at enrichment/depletion of orthogroup categories or interpro domains. 

#see other .R files for analyses related WGCNA, hierarchical clustering, and species-specific amino acid replacements.

#load libraries
library(airway)
library(tximport)
library(DESeq2)
library(PoiClaClu)
library(reshape)


#set working directory

setwd("/Users/gavin/genome/transcriptome_revisions_12-2023/all_genes/")
dir <- "/Users/gavin/genome/transcriptome_revisions_12-2023/all_genes/quants"

#get sample data in there
	#this is at https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/inopinata_sample_table.tsv

coldata <- read.table("inopinata_sample_table.tsv", sep="\t", header=T)
#tell R where the transcript count data are
	#these are in https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/tree/main/G3_revisions_1/all_genes_DESeq2/quants/inopinata
coldata$files <- file.path(dir,"inopinata", coldata$sample_id)
files <- file.path(dir, "inopinata", coldata$sample_id)

#get the tx2gene file in there that DESeq2 needs
	#this is here https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/tx2gene_inopinata.txt
tx2gene <- read.table("tx2gene_inopinata.txt", sep="\t", header=T)

#import the Salmon transcript count data
txi <- tximport(coldata$files, type = "salmon", tx2gene = tx2gene)

#assign stages to the transcript and sample data
txi$group <- factor(txi$Stage,levels=c("L3","L4","adult"))
coldata$group <- factor(coldata$Stage,levels=c("L3","L4","adult"))

#alright, create a DESeqDataSet object from the Salmon transcript count data
dat <- DESeqDataSetFromTximport(txi, coldata, design= ~ group)

#get a matrix of count data from the DESeqDataSet object
countdata <- round(assays(dat)[["counts"]])

#get the number of genes
nrow(dat)
#[1] 21472

#remove genes with 0 counts
keep <- rowSums(counts(dat)) > 1
dat <- dat[keep,]

#how many genes now?
nrow(dat)
#[1] 19669
	#there are some genes with no transcripts


#regularized log transformation of transcript counts
rld <- rlog(dat, blind = FALSE)
head(assay(rld,assaywithDimnames=TRUE), 3)

#                         1        2          3          4         5         6
#Sp34_10000100.t1 6.2878511 6.362370  1.7266236 3.39006838 5.2213254 5.3031134
#Sp34_10000300.t1 0.6653593 2.127672 -0.1559062 0.06796895 0.2100709 0.1416051
#Sp34_10000510.t1 1.2362822 1.413154  1.2509851 0.94206954 0.9463868 4.6584764
#                         7         8          9        10         11       12
#Sp34_10000100.t1 5.6876464 4.9971810  3.1005523 5.6012739  3.1317986 7.188805
#Sp34_10000300.t1 0.1436244 0.2531402 -0.6281303 0.4039753 -0.3031472 1.092456
#Sp34_10000510.t1 4.4744582 4.3973673  3.2348439 5.2982437  4.2876422 4.839365
#                       13        14
#Sp34_10000100.t1 4.103631 4.8921880
#Sp34_10000300.t1 1.041756 0.2437304
#Sp34_10000510.t1 5.302447 5.0069441

#write transformed data to file

write.table(assay(rld), file="inopinata_all_genes_rlog_transformed_transcript_counts.tsv", sep="\t", quote = FALSE)

#get factor levels right, default is alphabetical
dat$Stage <- factor(dat$Stage,levels=c("L3","L4","adult"))

#make a full model for getting at the relationship between transcript abundance and developmental stage
ddsTC <- DESeqDataSet(dat, ~ Date_RNA_prepared + Num_stage)

#do the LRT with the reduced model
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ + Date_RNA_prepared)
#get the results
resTC <- results(ddsTC)
#results as data frame
resTC_df <- as.data.frame(resTC)
#adjust p value with BH method
resTC_df$p_adj_bh <- p.adjust(resTC_df$pvalue,method="BH")
#"is the relationship between transcript abundance and developmental time significant?" THat is, is adjusted p < 0.05?
resTC_df$is_sig <- ifelse(resTC_df$p_adj_bh < 0.05, TRUE,FALSE)

#sort results by the slope, either positive or negative
resTC_increasing <- resTC_df[order(-resTC_df$log2FoldChange),]
resTC_decreasing <- resTC_df[order(resTC_df$log2FoldChange),]
#sort results by p value
resTC_by_p <- resTC_df[order(resTC_df$p_adj_bh),]
#add a category-- is slope positive or negative
resTC_df$pos_neg <- ifelse(resTC$log2FoldChange < 0, "negative","positive")
#add gene name column 
resTC_df$gene <- rownames(resTC_df)
#write linear model results to a file
write.table(resTC_df, file= 'inopinata_time_course_LRT_results.tsv',sep="\t", quote = FALSE)

#make some tables for the supplement

resTC_dat <- read.table("inopinata_time_course_LRT_results.tsv", sep="\t", header=TRUE)
#genes with positive relationships
inop_tx_dat_coef_positive <- resTC_dat[resTC_dat$pos_neg == "positive",]
#genes with negative relationships
inop_tx_dat_coef_negative <- resTC_dat[resTC_dat$pos_neg == "negative",]
#Supplemental table 8
write.table(inop_tx_dat_coef_negative[order(inop_tx_dat_coef_negative$p_adj_bh),], file= 'inopinata_time_course_LRT_results_negative.tsv',sep="\t", quote = FALSE)
#Supplemental table 9
write.table(inop_tx_dat_coef_positive[order(inop_tx_dat_coef_positive$p_adj_bh),], file= 'inopinata_time_course_LRT_results_positive.tsv',sep="\t", quote = FALSE)

#make some plots

#read transformed transcript data
#this was made above, line 73, also available here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/inopinata_all_genes_rlog_transformed_transcript_counts.tsv
inop_tx_dat <- read.table("inopinata_all_genes_rlog_transformed_transcript_counts.tsv", sep="\t", header=TRUE)

#sample data
	#here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/inopinata_sample_table.tsv
metadata <- read.table("inopinata_sample_table.tsv",header=TRUE,sep='\t')
#get factor levels right, default is alphabetical
metadata$Stage <- factor(metadata$Stage,levels = c("L3", "L4", "adult"))

#merge tx count data with linear model results 
inop_tx_dat_coef <- merge(inop_tx_dat,resTC_dat)
#sometimes, the above function does not work well, that data is also available here: 
	#https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/inop_tx_dat_coef.tsv
inop_tx_dat_coef <- read.table("inop_tx_dat_coef.tsv", sep="\t", header=TRUE)

#melt the data
meltdat <- melt(inop_tx_dat_coef,id.var=c("gene","is_sig","pos_neg","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","p_adj_bh"))

#replace columns with vague names
names(meltdat)[names(meltdat) == 'variable'] <- 'Sample_ID'
names(meltdat)[names(meltdat) == 'value'] <- 'expression'

#assign stages!
meltdat$stage[meltdat$Sample_ID == "s1"] <- "L3"
meltdat$stage[meltdat$Sample_ID == "s2"] <- "L3"
meltdat$stage[meltdat$Sample_ID == "s3"] <- "L3"
meltdat$stage[meltdat$Sample_ID == "s4"] <- "L3"
meltdat$stage[meltdat$Sample_ID == "s5"] <- "L3"
meltdat$stage[meltdat$Sample_ID == "s6"] <- "L4"
meltdat$stage[meltdat$Sample_ID == "s7"] <- "L4"
meltdat$stage[meltdat$Sample_ID == "s9"] <- "L4"
meltdat$stage[meltdat$Sample_ID == "s10"] <- "L4"
meltdat$stage[meltdat$Sample_ID == "s11"] <- "adult"
meltdat$stage[meltdat$Sample_ID == "s12"] <- "adult"
meltdat$stage[meltdat$Sample_ID == "s13"] <- "adult"
meltdat$stage[meltdat$Sample_ID == "s14"] <- "adult"
meltdat$stage[meltdat$Sample_ID == "s15"] <- "adult"
#get a "no_change" value in there
meltdat$pos_neg[meltdat$is_sig == FALSE] <- "no_change"


#next, for plotting, the aim is to get the average abundance for each gene at each stage. we're also going to plot genes with postive, negative, or no obvious relationship with time in different facets. so, we're going to do this.

#make a column with the stage, gene, and pos/neg info
meltdat$stage.gene.is_sig.pos_neg <- paste(meltdat$stage,meltdat$gene,meltdat$is_sig,meltdat$pos_neg)

#because each gene-stage combination has the same value for this "stage.gene.is_sig.pos_neg" column, this next line will find the average across all samples with the same value of "stage.gene.is_sig.pos_neg" -- this will give us the numbers we need (average expression for each gene at each stage) for the plot we want.
meltdat_merge_agg <- aggregate(expression ~ stage.gene.is_sig.pos_neg, FUN=mean,data=meltdat)

#load a package we need to _split_ the combined column we just made to get averages

library(tidyr)

#split into two
meltdat_merge_agg<- separate(meltdat_merge_agg, stage.gene.is_sig.pos_neg, into = c("stage.gene.is_sig","pos_neg"), sep = " (?=[^ ]+$)")

#split into two
meltdat_merge_agg<- separate(meltdat_merge_agg, stage.gene.is_sig, into = c("stage.gene","is_sig"), sep = " (?=[^ ]+$)")

#split into two-- we're done splitting columns
meltdat_merge_agg<- separate(meltdat_merge_agg, stage.gene, into = c("stage","gene"), sep = " (?=[^ ]+$)")

#check to see if this worked-- should be one row for each stage. "expression" is mean expression among samples at that stage for that gene
meltdat_merge_agg[meltdat_merge_agg$gene=="Sp34_10000100.t1",]
	#yay it worked

#okay, load some packages for making plots
library(ggplot2)
library(cowplot)
library(lemon)


#get factor levels right for plotting
meltdat_merge_agg$stage <- factor(meltdat_merge_agg$stage,levels = c("L3", "L4", "adult"))

meltdat_merge_agg$pos_neg[meltdat_merge_agg$pos_neg == 'negative'] <- 'Negative'
meltdat_merge_agg$pos_neg[meltdat_merge_agg$pos_neg == 'no_change'] <- 'No change'
meltdat_merge_agg$pos_neg[meltdat_merge_agg$pos_neg == 'positive'] <- 'Positive'



panel_a <- ggplot(meltdat_merge_agg, aes(x=stage, y=expression,group=gene)) + geom_line(alpha=0.05,size=0.05,aes(colour=is_sig))  + facet_rep_wrap(~pos_neg) + scale_color_manual(values=c("black","red")) + theme_cowplot() + ylab("Transformed expression") + xlab("Stage") + theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.015,size=24),strip.background=element_blank()) + ggtitle("a")
	#this is supplemental figure 3A




#ORTHOLOG CATEGORIES

#original orthofinder run code on lines 242-343 of pipeline.sh, here: (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/pipeline.sh)

#then, those orthologous groups were classified for the revision. that code is here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/orthologous_genes_notes.sh

#and here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/join_orthogroups_genes_categories_notes.R

#the file inopinata_GENES_ogs_categories.tsv is here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/inopinata_GENES_ogs_categories.tsv

#get the ortholog category-gene id data in there
inop_ortho <- read.table("inopinata_GENES_ogs_categories.tsv", sep="\t", header=FALSE)
#rename column
inop_ortho$gene <- inop_ortho$V1
#combine transcript model coefficients per gene with ortholog information
ortho_merge <- merge(inop_tx_dat_coef,inop_ortho, by.x="gene", by.y="gene")
#get a "no_change label"
ortho_merge$pos_neg[ortho_merge$is_sig == FALSE] <- "no_change"

#get the counts of ortholog categories for ALL genes in the genome
total_ortho_cat <- as.data.frame(table(inop_ortho$V3))
#give it a label
total_ortho_cat$pos_neg <- "All genes"

#separate df's for each transcript dynamics category
pos <- ortho_merge[ortho_merge$pos_neg == "positive",]
neg <- ortho_merge[ortho_merge$pos_neg == "negative",]
noc <- ortho_merge[ortho_merge$pos_neg == "no_change",]

#get the nubmers of the ortholog categories (single-copy, inopinata orphan, etc.) for each transcript dyanmics category (positive, negative, no change)
pos_cat <- as.data.frame(table(pos$V3))
neg_cat <- as.data.frame(table(neg$V3))
noc_cat <- as.data.frame(table(noc$V3))

#give them a good label for plotting
pos_cat$pos_neg <- "Positive"
neg_cat$pos_neg <- "Negative"
noc_cat$pos_neg <- "No change"


#add all of the dfs for each tx category into one big df for plotting
ortho_cat_df <- rbind(total_ortho_cat,pos_cat,neg_cat,noc_cat)

#get the label factor right
ortho_cat_df$pos_neg <- as.factor(ortho_cat_df$pos_neg)

ortho_cat_df$pos_neg <- factor(ortho_cat_df$pos_neg,levels = c("All genes","Negative","No change","Positive"))

#get the total number of genes for each category to print the sample sizes on the plot
tot_genes <- aggregate(Freq ~ pos_neg, FUN=sum, data=ortho_cat_df)

#load a useful library for making composite plots
library(patchwork)

#this is supplemental figure 3B
panel_b <- ggplot(ortho_cat_df, aes(x = pos_neg, y = Freq, fill = Var1)) + geom_col(position = "fill") + scale_fill_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="bottom",legend.title=element_blank(),legend.justification = "center",plot.title = element_text(face = "bold",hjust = -0.015,size=24)) + annotate("text",x=1, y=0.05, label=tot_genes$Freq[1],size=5) + annotate("text",x=2, y=0.05, label=tot_genes$Freq[2],size=5) + annotate("text",x=3, y=0.05, label=tot_genes$Freq[3],size=5) + annotate("text",x=4, y=0.05, label=tot_genes$Freq[4],size=5) + ylab("Proportion") + ggtitle("b") + xlab("Transcriptional dynamics pattern") + guides(fill=guide_legend(nrow=2,byrow=TRUE))

#this is supplemental figure 3B YEAH are you ready for this reviewers and anyone else who cares about this
panel_a/panel_b
	#okay, yes, I cleaned up the figure legend in adobe illustrator

#saved it as some files
#ggsave("all_inopinata_genes_by_developmental_dynamics_line_stacked_bar_iii.png",units="in",height=7,width=10)
#ggsave("all_inopinata_genes_by_developmental_dynamics_line_stacked_bar_iii.pdf",units="in",height=7,width=10,useDingbats=TRUE)




#okay, chi-square (chi square) for comparing the proportion of OG categories in each of our transcriptional dynamic categories (positive, negative, no change) with those of the whole protein set

#load this library if we haven't already
library(reshape2)

#rename factor level 
levels(ortho_cat_df$pos_neg)[match("All genes",levels(ortho_cat_df$pos_neg))] <- "All_genes"
#CAST-- the idea is to make a data frame for OG category chi square
orthocast <- dcast(ortho_cat_df,Var1~pos_neg,value.var="Freq")

#rename some columns
names(orthocast)[names(orthocast) == 'Negative'] <- 'negative'
names(orthocast)[names(orthocast) == 'No change'] <- 'no_change'
names(orthocast)[names(orthocast) == 'Positive'] <- 'positive'

#get all genes sums again
tot_genes <- aggregate(Freq ~ pos_neg, FUN=sum, data=ortho_cat_df)

#how many genes _aren't_ in some category. this requires the tot_genes df above
orthocast$all_genes_no <- tot_genes$Freq[1]-orthocast$All_genes
orthocast$positive_no <- tot_genes$Freq[4]-orthocast$positive
orthocast$no_change_no <- tot_genes$Freq[3]-orthocast$no_change
orthocast$negative_no <- tot_genes$Freq[2]-orthocast$negative

#subset just the genes related to positive trajectories over time
poscast <- data.frame(orthotype=orthocast$Var1, All_genes=orthocast$All_genes,All_genes_no= orthocast$all_genes_no,positive=orthocast$positive,positive_no=orthocast$positive_no)

#make an empty variable that we'll add stuff (ie chi square stats) to for the loop
posdf <- NULL
	#for rows in dataframe do
for(i in 1:nrow(poscast)){
	#get a row
	orthoclass <- poscast$orthotype[i]
	#make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that are and are not that specific OG category among positive trajectory genes and all genes in the genome
	the_matrix <- matrix(c(poscast[i,2],poscast[i,3],poscast[i,4],poscast[i,5]),ncol=2)
	#do a chi square test on that matrix
	matrixtest <- chisq.test(the_matrix)
	#make a row for our df including:
		#the orthogroup category
		#the number of genes in the genome that are in that category
		#the number of genes in the genome that are not in that category
		#the number of positive trajectory genes in that category
		#the number of positive trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of genes among those with a positive trajectory in that category
		#the chi square test statistic
		#the chisquare p value
	row_to_add <- data.frame(orthotype=orthoclass,all_genes_with=poscast[i,2],all_genes_without=poscast[i,3], test_genes_with=poscast[i,4],test_genes_without=poscast[i,5], fra_all_genes_with=poscast[i,2]/21442, fra_test_genes_with=poscast[i,4]/6658, chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value)
	#add that row to the df and repeat
	posdf <- rbind(posdf,row_to_add)
}

#empty variables for "no change" and "negative" trajectories
nochangedf <- NULL
negdf <- NULL

#data frames for no change and negative trajectories
nochangecast <- data.frame(orthotype=orthocast$Var1, All_genes=orthocast$All_genes,All_genes_no= orthocast$all_genes_no,no_change=orthocast$no_change,no_change_no=orthocast$no_change_no)
negcast <- data.frame(orthotype=orthocast$Var1, All_genes=orthocast$All_genes,All_genes_no= orthocast$all_genes_no,negative=orthocast$negative,negative_no=orthocast$negative_no)

#continuing chi square stats
#no change,
		#for rows in dataframe do
for(i in 1:nrow(nochangecast)){
	#get a row
	orthoclass <- nochangecast$orthotype[i]
	#make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that are and are not that specific OG category among no change trajectory genes and all genes in the genome
	the_matrix <- matrix(c(nochangecast[i,2],nochangecast[i,3],nochangecast[i,4],nochangecast[i,5]),ncol=2)
	#do a chi square test on that matrix
	matrixtest <- chisq.test(the_matrix)
	#make a row for our df including:
		#the orthogroup category
		#the number of genes in the genome that are in that category
		#the number of genes in the genome that are not in that category
		#the number of no change trajectory genes in that category
		#the number of no change trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of genes among those with a no change trajectory in that category
		#the chi square test statistic
		#the chisquare p value
	row_to_add <- data.frame(orthotype=orthoclass,all_genes_with=nochangecast[i,2],all_genes_without=nochangecast[i,3], test_genes_with=nochangecast[i,4],test_genes_without=nochangecast[i,5], fra_all_genes_with=nochangecast[i,2]/21442, fra_test_genes_with=nochangecast[i,4]/8427, chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value)
	#add that row to the df and repeat
	nochangedf <- rbind(nochangedf,row_to_add)
}

#continuing chi square stats
#negative trajectories,
		#for rows in dataframe do

for(i in 1:nrow(negcast)){
	#get a row
	orthoclass <- negcast$orthotype[i]
	#make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that are and are not that specific OG category among negative trajectory genes and all genes in the genome
	the_matrix <- matrix(c(negcast[i,2],negcast[i,3],negcast[i,4],negcast[i,5]),ncol=2)
	#do a chi square test on that matrix
	matrixtest <- chisq.test(the_matrix)
	#make a row for our df including:
		#the orthogroup category
		#the number of genes in the genome that are in that category
		#the number of genes in the genome that are not in that category
		#the number of negative trajectory genes in that category
		#the number of negative trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of genes among those with a negative trajectory in that category
		#the chi square test statistic
		#the chisquare p value	
		row_to_add <- data.frame(orthotype=orthoclass,all_genes_with=negcast[i,2],all_genes_without=negcast[i,3], test_genes_with=negcast[i,4],test_genes_without=negcast[i,5], fra_all_genes_with=negcast[i,2]/21442, fra_test_genes_with=negcast[i,4]/4426, chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value)
		#add that row to the df and repeat
	negdf <- rbind(negdf,row_to_add)
}

#add columns to note tx category
posdf$tx_direction <- "positive"
nochangedf$tx_direction <- "no change"
negdf$tx_direction <- "negative"
#put the ortholog category chi square tests in a single df
all_orthoclass_chisquare <- rbind(posdf,nochangedf,negdf)
#adjust p values for multiple testing
all_orthoclass_chisquare$p.adjust <- p.adjust(all_orthoclass_chisquare$chisq_p,method="BH")

#write this to a file-- this is a supplemental table! supplemental table 12 !!
write.table(all_orthoclass_chisquare,"all_genes_ortholog_type_tx_direction_chi_square.tsv",quote=FALSE,row.names = FALSE,sep="\t")

#next, domain info

#ran interproscan on protein sets to get domain info. that code is here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/interproscan/workflow_notes.sh
#the file inopinata_04.tsv is here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/interproscan/inopinata_04.tsv


#get inopinata domain info

inopinata_interpro <- read.table("inopinata_04.tsv", sep="\t")

#make column names
names(inopinata_interpro)[names(inopinata_interpro) == 'V1'] <- 'gene'
names(inopinata_interpro)[names(inopinata_interpro) == 'V2'] <- 'interpro_id'
names(inopinata_interpro)[names(inopinata_interpro) == 'V3'] <- 'interpro_description'
#merge with expression data
big_df_merge <- merge(inopinata_interpro,ortho_merge)

#count the domains (this is the number of genes that carry each domain)

inopinata_domain_counts <- as.data.frame(table(inopinata_interpro$interpro_id))

#get just id's and descriptions

interpro_domains_descriptions <- data.frame(interpro_id=inopinata_interpro$interpro_id,interpro_description=inopinata_interpro$interpro_description)

#ensure there are no duplicate entries

interpro_domains_descriptions <- interpro_domains_descriptions[!duplicated(interpro_domains_descriptions), ]

#make readable column names
names(inopinata_domain_counts)[names(inopinata_domain_counts) == 'Var1'] <- 'interpro_id'
names(inopinata_domain_counts)[names(inopinata_domain_counts) == 'Freq'] <- 'count'

#merge inopinata domain counts with interpro domain descriptions
inopinata_domain_counts <- merge(inopinata_domain_counts,interpro_domains_descriptions)

#get the sum of the counts
inopinata_domain_counts$total <- sum(inopinata_domain_counts$count)
#get the number of genes with domains that ARE NOT that domain
inopinata_domain_counts$num_not_domain <- inopinata_domain_counts$total-inopinata_domain_counts$count
#get the fraction
inopinata_domain_counts$fra_domain <- inopinata_domain_counts$count/inopinata_domain_counts$total
	#here, I believe the total is the total number of unique gene:domain relationships in the genome. Not allowing duplicate domains _per gene_.
sum(inopinata_domain_counts$count)
#[1] 31419
#this is the total number of unique gene:domain relationships

#add label for comparing the total gene set with the transcriptional categories (positive, negative, no change)
inopinata_domain_counts$pos_neg <- "All genes"

#okay, count the number of genes in big_df_merge (these are the genes with gene:domain relationships)
length(unique(big_df_merge$gene))
#[1] 8963

#get the different tx categories
pos <- big_df_merge[big_df_merge$pos_neg == "positive",]
neg <- big_df_merge[big_df_merge$pos_neg == "negative",]
noc <- big_df_merge[big_df_merge$pos_neg == "no_change",]
#get domain counts
pos_cat <- as.data.frame(table(pos$interpro_id))
neg_cat <- as.data.frame(table(neg$interpro_id))
noc_cat <- as.data.frame(table(noc$interpro_id))
#add tx categories to domain counts
pos_cat$pos_neg <- "positive"
neg_cat$pos_neg <- "negative"
noc_cat$pos_neg <- "no_change"

#rename column ids
names(pos_cat)[names(pos_cat) == 'Var1'] <- 'interpro_id'
names(pos_cat)[names(pos_cat) == 'Freq'] <- 'count'

names(neg_cat)[names(neg_cat) == 'Var1'] <- 'interpro_id'
names(neg_cat)[names(neg_cat) == 'Freq'] <- 'count'

names(noc_cat)[names(noc_cat) == 'Var1'] <- 'interpro_id'
names(noc_cat)[names(noc_cat) == 'Freq'] <- 'count'


#okay, chi-square (chi square) for comparing the proportion of interpro domains in each of our transcriptional dynamic categories (positive, negative, no change) with those of the whole protein set

#positive trajectory genes

#select columns I want (interpro id, domain count, trajectory id)
inop_dom_count_abb <- data.frame(interpro_id=inopinata_domain_counts$interpro_id,count=inopinata_domain_counts$count,pos_neg=inopinata_domain_counts$pos_neg)
#get total number of gene:domain pairs for whole genome
inop_dom_count_abb$total_all_genes <- sum(inop_dom_count_abb$count)
#get difference between total number of gene:domain pairs and the value for each specific domain for whole genome
inop_dom_count_abb$num_all_genes_not_domain <- inop_dom_count_abb$total_all_genes-inop_dom_count_abb$count

#get total number of gene:domain pairs for positive tx genes
pos_cat$total_test_genes <- sum(pos_cat$count)
#get difference between total number of gene:domain pairs and the value for each specific domain for positive tx genes
pos_cat$num_test_genes_not_domain <- pos_cat$total_test_genes-pos_cat$count
#rename for unique column id
names(pos_cat)[names(pos_cat) == 'count'] <- 'pos_count'

#load library if unloaded
library(dplyr)
#remove unnecessary column
pos_cat <- pos_cat %>% select(-pos_neg)

#rename for unique column id
names(inop_dom_count_abb)[names(inop_dom_count_abb) == 'count'] <- 'all_count'
#remove unnecessary column
inop_dom_count_abb <- inop_dom_count_abb %>% select(-pos_neg)

#merge positive tx trajectory and whole genome df's
pos_tot_merge <- merge(inop_dom_count_abb,pos_cat, all.x=TRUE)

#replace NA's with appropriate values (i.e., zero or total number of pos genes)
pos_tot_merge$pos_count[is.na(pos_tot_merge$pos_count)] <- 0
pos_tot_merge$total_test_genes[is.na(pos_tot_merge$total_test_genes)] <- sum(pos_cat$pos_count)
#refigure num_test_genes_not_domain
pos_tot_merge$num_test_genes_not_domain <- pos_tot_merge$total_test_genes-pos_tot_merge$pos_count


#okay, let's try some chi square (as before with OG categories)
poscast <- data.frame(interpro_id=pos_tot_merge$interpro_id, All_genes=pos_tot_merge$all_count,All_genes_no= pos_tot_merge$num_all_genes_not_domain,positive=pos_tot_merge$pos_count,positive_no=pos_tot_merge$num_test_genes_not_domain)

#make an empty variable that we'll add stuff (ie chi square stats) to for the loop
posdf <- NULL
	#for rows in dataframe do
for(i in 1:nrow(poscast)){
	#get a row
	interpro_id <- poscast$interpro_id[i]
	#make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that are and are not that specific interpro domain among positive trajectory genes and all genes in the genome
	the_matrix <- matrix(c(poscast[i,2],poscast[i,3],poscast[i,4],poscast[i,5]),ncol=2)
	#do a chi square test on that matrix
	matrixtest <- chisq.test(the_matrix)
	#make a row for our df including:
		#the orthogroup category
		#the number of genes in the genome that are in that category
		#the number of genes in the genome that are not in that category
		#the number of positive trajectory genes in that category
		#the number of positive trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of genes among those with a positive trajectory in that category
		#the chi square test statistic
		#the chisquare p value
	row_to_add <- data.frame(interpro_id=interpro_id,all_genes_with=poscast[i,2],all_genes_without=poscast[i,3], test_genes_with=poscast[i,4],test_genes_without=poscast[i,5], fra_all_genes_with=poscast[i,2]/31419, fra_test_genes_with=poscast[i,4]/10508, chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value)
	#add that row to the df and repeat
	posdf <- rbind(posdf,row_to_add)
}


#negative trajectory genes chi square chi-square


#get total number of gene:domain pairs for negative tx genes
neg_cat$total_test_genes <- sum(neg_cat$count)
#get difference between total number of gene:domain pairs and the value for each specific domain for negative tx genes
neg_cat$num_test_genes_not_domain <- neg_cat$total_test_genes-neg_cat$count
#rename for unique column id
names(neg_cat)[names(neg_cat) == 'count'] <- 'neg_count'

#remove unnecessary column
neg_cat <- neg_cat %>% select(-pos_neg)

#merge negative and whole genome df's
neg_tot_merge <- merge(inop_dom_count_abb,neg_cat, all.x=TRUE)

#replace NA's with appropriate values (i.e., zero or total number of pos genes)
neg_tot_merge$neg_count[is.na(neg_tot_merge$neg_count)] <- 0
neg_tot_merge$total_test_genes[is.na(neg_tot_merge$total_test_genes)] <- sum(neg_cat$neg_count)
#refigure num_test_genes_not_domain
neg_tot_merge$num_test_genes_not_domain <- neg_tot_merge$total_test_genes-neg_tot_merge$neg_count

#okay, let's try some chi square
poscast <- data.frame(interpro_id=neg_tot_merge$interpro_id, All_genes=neg_tot_merge$all_count,All_genes_no= neg_tot_merge$num_all_genes_not_domain,negative=neg_tot_merge$neg_count,negative_no=neg_tot_merge$num_test_genes_not_domain)


#make an empty variable that we'll add stuff (ie chi square stats) to for the loop
negdf <- NULL
	#for rows in dataframe do
for(i in 1:nrow(poscast)){
	#get a row
	interpro_id <- poscast$interpro_id[i]
	#make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that are and are not that specific interpro domain among negative trajectory genes and all genes in the genome
	the_matrix <- matrix(c(poscast[i,2],poscast[i,3],poscast[i,4],poscast[i,5]),ncol=2)
	#do a chi square test on that matrix
	matrixtest <- chisq.test(the_matrix)
	#make a row for our df including:
		#the orthogroup category
		#the number of genes in the genome that are in that category
		#the number of genes in the genome that are not in that category
		#the number of negative trajectory genes in that category
		#the number of negative trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of genes among those with a negative trajectory in that category
		#the chi square test statistic
		#the chisquare p value
	row_to_add <- data.frame(interpro_id=interpro_id,all_genes_with=poscast[i,2],all_genes_without=poscast[i,3], test_genes_with=poscast[i,4],test_genes_without=poscast[i,5], fra_all_genes_with=poscast[i,2]/31419, fra_test_genes_with=poscast[i,4]/7827, chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value)
	#add that row to the df and repeat
	negdf <- rbind(negdf,row_to_add)
}

#no change genes chi square


#get total number of gene:domain pairs for no change tx genes
noc_cat$total_test_genes <- sum(noc_cat$count)
#get difference between total number of gene:domain pairs and the value for each specific domain for no change tx genes
noc_cat$num_test_genes_not_domain <- noc_cat$total_test_genes-noc_cat$count
#rename for unique column id
names(noc_cat)[names(noc_cat) == 'count'] <- 'noc_count'

#remove unnecessary column
noc_cat <- noc_cat %>% select(-pos_neg)

#merge positive and whole genome df's
noc_tot_merge <- merge(inop_dom_count_abb,noc_cat, all.x=TRUE)

#replace NA's with appropriate values (i.e., zero or total number of pos genes)
noc_tot_merge$noc_count[is.na(noc_tot_merge$noc_count)] <- 0
noc_tot_merge$total_test_genes[is.na(noc_tot_merge$total_test_genes)] <- sum(noc_cat$noc_count)
#refigure num_test_genes_not_domain
noc_tot_merge$num_test_genes_not_domain <- noc_tot_merge$total_test_genes-noc_tot_merge$noc_count

#okay, let's try some chi square
poscast <- data.frame(interpro_id=noc_tot_merge$interpro_id, All_genes=noc_tot_merge$all_count,All_genes_no= noc_tot_merge$num_all_genes_not_domain,noc=noc_tot_merge$noc_count,noc_no=noc_tot_merge$num_test_genes_not_domain)

#make an empty variable that we'll add stuff (ie chi square stats) to for the loop
nocdf <- NULL
	#for rows in dataframe do
for(i in 1:nrow(poscast)){
	#get a row
	interpro_id <- poscast$interpro_id[i]
	#make a matrix from that row, specifically, a 2x2 contingency table getting at the genes that are and are not that specific interpro domain among no change trajectory genes and all genes in the genome
	the_matrix <- matrix(c(poscast[i,2],poscast[i,3],poscast[i,4],poscast[i,5]),ncol=2)
	#do a chi square test on that matrix
	matrixtest <- chisq.test(the_matrix)
	#make a row for our df including:
		#the orthogroup category
		#the number of genes in the genome that are in that category
		#the number of genes in the genome that are not in that category
		#the number of no change trajectory genes in that category
		#the number of no change trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of genes among those with a no change trajectory in that category
		#the chi square test statistic
		#the chisquare p value
	row_to_add <- data.frame(interpro_id=interpro_id,all_genes_with=poscast[i,2],all_genes_without=poscast[i,3], test_genes_with=poscast[i,4],test_genes_without=poscast[i,5], fra_all_genes_with=poscast[i,2]/31419, fra_test_genes_with=poscast[i,4]/11654, chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value)
	#add that row to the df and repeat
	nocdf <- rbind(nocdf,row_to_add)
}



#join to the interpro descriptions

posdf <- merge(posdf, interpro_domains_descriptions)
negdf <- merge(negdf, interpro_domains_descriptions)
nocdf <- merge(nocdf, interpro_domains_descriptions)


#add labels back, so we can adjust for multiple testing, then split again
posdf$pos_neg <- "positive"
negdf$pos_neg <- "negative"
nocdf$pos_neg <- "no_change"

all_chisq_df <- rbind(posdf,negdf,nocdf)

#adjust p for multiple testing
all_chisq_df$p.adj <- p.adjust(all_chisq_df$chisq_p, method="BH")

#split up again
posdf <- all_chisq_df[all_chisq_df$pos_neg == "positive",]
negdf <- all_chisq_df[all_chisq_df$pos_neg == "negative",]
nocdf <- all_chisq_df[all_chisq_df$pos_neg == "no_change",]


#for some reason "write.table" is not doing what i want?
	#this helped
posdf <- posdf %>% mutate_if(is.character, trimws)
negdf <- negdf %>% mutate_if(is.character, trimws)
nocdf <- nocdf %>% mutate_if(is.character, trimws)

#this is supplemental table 10 (positive trajectory genes)
write.table(posdf,"all_genes_DOMAIN_tx_positive_chi_square_ii.tsv",quote=FALSE,row.names = FALSE,sep="\t")

#this is supplemental table 11 (negative trajectory genes)
write.table(negdf,"all_genes_DOMAIN_tx_negative_chi_square_ii.tsv",quote=FALSE,row.names = FALSE,sep="\t")

#not in the supplement but is included on github (no change genes)
write.table(nocdf,"all_genes_DOMAIN_tx_no_change_chi_square_ii.tsv",quote=FALSE,row.names = FALSE,sep="\t")


















