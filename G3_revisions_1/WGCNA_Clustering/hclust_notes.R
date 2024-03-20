##This is the code for the hierarchical clustering analyses, with all C. inopinata genes (not just single-copy orthologs), for Woodruff et al. "Widespread changes in gene expression accompany body size evolution in nematodes"

#load libraries

library(reshape2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(lemon)
library(dplyr)

#set working directory
setwd("/Users/gavin/genome/transcriptome_revisions_12-2023/WGCNA/")

#get data in there...

	#transformed expression data; generated in DESeq2_all_genes_notes.R , https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/DESeq2_all_genes_notes.R

inop_tx_dat  <- read.table("inopinata_all_genes_rlog_transformed_transcript_counts.tsv", sep="\t", header=TRUE,row.names=1)
#gene column
inop_tx_dat$gene <- rownames(inop_tx_dat)

#get transcript data in there in matrix form for clustering
inop_tx_dat_mat <- as.matrix(read.table("inopinata_all_genes_rlog_transformed_transcript_counts.tsv", sep="\t", header=TRUE,row.names=1))

#transpose
tdat = t(inop_tx_dat_mat)

#get distances
ino_dat_distm <- dist(inop_tx_dat_mat)
#do hierarchical clustering
gene_hclust <- hclust(ino_dat_distm, method = "complete")


#elbow method for picking some number of clusters
#mydata <- inop_tx_dat_mat
#wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  #for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",
     #ylab="Within groups sum of squares")
#
#
#plot(gene_hclust, labels=FALSE)

#4-10 clusters looks reasonable by elbow method, but I chose k=20 to make it at least somewhat more comparable to the WGCNA analysis (which had 52 modules)


#let's try 20 clusters, slice the HC tree into 20 clusters
clust_k_20 <- cutree(gene_hclust, k = 20)

#get factors right
inop_tx_dat$clust_k_20 <- as.factor(clust_k_20)


#melt data
meltdat <- melt(inop_tx_dat,id.var=c("gene","clust_k_20"))
#replace column ids
names(meltdat)[names(meltdat) == 'variable'] <- 'Sample_ID'
names(meltdat)[names(meltdat) == 'value'] <- 'expression'
#get sample data in there;  inopinata_sample_table.tsv, can be found here: https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/all_genes_DESeq2/inopinata_sample_table.tsv

metadata <- read.table("inopinata_sample_table.tsv",header=TRUE,sep='\t')
#get stage factor levels right
metadata$Stage <- factor(metadata$Stage,levels = c("L3", "L4", "adult"))
#merge transcript data with sample metadata
meltdat_merge <- merge(meltdat,metadata)


#get a stage-gene-cluster column for getting the means of these categories for making the plots we want
meltdat_merge$stage.gene.clust <- paste(meltdat_merge$Stage,meltdat_merge$gene,meltdat_merge$clust_k_20)
#get averages for each gene at each given stage
meltdat_merge_agg <- aggregate(expression ~ stage.gene.clust, FUN=mean,data=meltdat_merge)
#break up the "stage.gene.clust" column into three columns
meltdat_merge_agg<- separate(meltdat_merge_agg, stage.gene.clust, into = c("stage.gene","cluster"), sep = " (?=[^ ]+$)")

meltdat_merge_agg<- separate(meltdat_merge_agg, stage.gene, into = c("stage","gene"), sep = " (?=[^ ]+$)")


#get factor levels right
meltdat_merge_agg$stage <- factor(meltdat_merge_agg$stage,levels = c("L3", "L4", "adult"))

meltdat_merge_agg$cluster <- factor(meltdat_merge_agg$cluster,levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"))

#make that clusters plot
ggplot(meltdat_merge_agg, aes(x=stage, y=expression,group=gene)) + geom_line(alpha=0.15,size=0.15) + facet_rep_wrap(~cluster) + scale_color_brewer(palette="Set1") + theme_cowplot() +ylab("Transformed expression (rlog)") + xlab("Stage")

	#this is supplemental figure 6

#save the clusters
gene_clust_key <- data.frame(gene=inop_tx_dat$gene,cluster=inop_tx_dat$clust_k_20)
write.table(gene_clust_key, "inopinata_genes_h_clusters.tsv",quote=FALSE,sep="\t")
	#this is supplemental table 17


#linear models for each cluster
#make an empty variable that we'll add stuff (model stats) to for the loop
mod_coeff <- NULL
#get factors right for clusters
meltdat_merge_agg$cluster <- as.factor(meltdat_merge_agg$cluster)

#for each cluster...
for (i in levels(meltdat_merge_agg$cluster)){
        #get genes just in that cluster
    guh <- meltdat_merge_agg[meltdat_merge_agg$cluster == i,]
        #get a data frame of OLS linear model coefficients
    buh <- as.data.frame(summary(lm(expression ~ stage,data=guh))$coefficients)
        #print the model coefficients to a row
    duh <- data.frame(L4_coeff = buh[2,1], adult_coeff = buh[3,1], L4_p=buh[2,4], adult_p=buh[3,4],cluster=i)
        #add row to growing data frame
    mod_coeff <- rbind(mod_coeff, duh)
} 


#adjust p values for multiple testing
mod_coeff$adult_p_adjust <- p.adjust(mod_coeff$adult_p,method="BH")

mod_coeff$L4_p_adjust <- p.adjust(mod_coeff$L4_p,method="BH")
#order by adult coefficients
mod_coeff[order(-mod_coeff$adult_coeff,mod_coeff$adult_p_adjust),]


#write to table
write.table(mod_coeff[order(-mod_coeff$adult_coeff,mod_coeff$adult_p_adjust),], "h_cluster_linear_model_coefficients.tsv",quote=FALSE,sep="\t")
	#this is supplemental table 18



#ok. tie to orthogroup categories for chi square

#the clusters with the most interesting dynamics are 20,15,8 (increasing); 19 (max L4); and 16, 17 (decreasing)
#get ortholog type categories
		#inopinata gene ortholog categories; generated in join_orthogroups_genes_categories_notes.R , https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/orthologs/join_orthogroups_genes_categories_notes.R


inop_ortho <- read.table("inopinata_GENES_ogs_categories.tsv", sep="\t", header=FALSE)
#get a good gene name column
inop_ortho$gene <- inop_ortho$V1
#merge transcript data with ortholog type data
ortho_merge <- merge(inop_tx_dat,inop_ortho, by.x="gene", by.y="gene")
#get ortholog type counts for all inopinata genes
total_ortho_cat <- as.data.frame(table(inop_ortho$V3))
total_ortho_cat$cluster <- "All genes"
#get the clusters of interest
	#the clusters with the most interesting dynamics are 20,15,8 (increasing); 19 (max L4); and 16, 17 (decreasing)
cl20 <- ortho_merge[ortho_merge$clust_k_20 == "20",]
cl15 <- ortho_merge[ortho_merge$clust_k_20 == "15",]
cl8 <- ortho_merge[ortho_merge$clust_k_20 == "8",]
cl19 <- ortho_merge[ortho_merge$clust_k_20 == "19",]
cl16 <- ortho_merge[ortho_merge$clust_k_20 == "16",]
cl17 <- ortho_merge[ortho_merge$clust_k_20 == "17",]
#count ortholog types
cl20_ortho_cat <- as.data.frame(table(cl20$V3))
cl15_ortho_cat <- as.data.frame(table(cl15$V3))
cl8_ortho_cat <- as.data.frame(table(cl8$V3))
cl19_ortho_cat <- as.data.frame(table(cl19$V3))
cl16_ortho_cat <- as.data.frame(table(cl16$V3))
cl17_ortho_cat <- as.data.frame(table(cl17$V3))
#put the cluster id in the df
cl20_ortho_cat$cluster <- "20"
cl15_ortho_cat$cluster <- "15"
cl8_ortho_cat$cluster <- "8"
cl19_ortho_cat$cluster <- "19"
cl16_ortho_cat$cluster <- "16"
cl17_ortho_cat$cluster <- "17"
#combine the data frames for plotting
ortho_cat_df <- rbind(total_ortho_cat,cl20_ortho_cat,cl15_ortho_cat,cl8_ortho_cat,cl19_ortho_cat,cl16_ortho_cat,cl17_ortho_cat)

#get factor right
ortho_cat_df$cluster <- as.factor(ortho_cat_df$cluster)

ortho_cat_df$cluster <- factor(ortho_cat_df$cluster,levels = c("All genes","8","15","20","19","16","17"))

#get total number of genes per cluster
tot_genes <- aggregate(Freq ~ cluster, FUN=sum, data=ortho_cat_df)


panel_b <- ggplot(ortho_cat_df, aes(x = cluster, y = Freq, fill = Var1)) + geom_col(position = "fill") + scale_fill_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="bottom",legend.title=element_blank(),legend.justification = "center",plot.title = element_text(face = "bold",hjust = -0.015,size=24)) + annotate("text",x=1, y=0.05, label=tot_genes$Freq[1],size=5) + annotate("text",x=2, y=0.05, label=tot_genes$Freq[2],size=5) + annotate("text",x=3, y=0.05, label=tot_genes$Freq[3],size=5) + annotate("text",x=4, y=0.05, label=tot_genes$Freq[4],size=5) + annotate("text",x=5, y=0.05, label=tot_genes$Freq[5],size=5) + annotate("text",x=6, y=0.05, label=tot_genes$Freq[6],size=5) + annotate("text",x=7, y=0.05, label=tot_genes$Freq[7],size=5) + xlab("Cluster") + ylab("Proportion") + ggtitle("b") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
	#this is supplemental figure 7B


#okay, chi square (chi-square) stats for the ortholog types

#replace "All genes" with "All_genes"
levels(ortho_cat_df$cluster)[match("All genes",levels(ortho_cat_df$cluster))] <- "All_genes"

#CAST that df to get a form I want for hypothesis tests
orthocast <- dcast(ortho_cat_df,Var1~cluster,value.var="Freq")
#replace NA with zero
orthocast[is.na(orthocast)] <- 0
#get empty variable

testfun<-NULL

#for all of the clusters,

for(i in levels(ortho_cat_df$cluster)[-1]){
#this is the cluster of interest
  the_i <- i
   #make an empty variable that we'll add stuff (ie chi square stats) the inner loop
  i_df <-NULL
  #make the df we need for doing chi square tests...
	#the ortholog category
	#the number of genes in the genome in those categories
	#the number of genes in the genome NOT in those categories
	#the number of genes in the cluster of interest in those categories
	#the number of genes in the cluster of interest NOT in those categories
  i_df <- data.frame(orthotype=orthocast$Var1, All_genes=orthocast[,2],All_genes_no= tot_genes[tot_genes$cluster == "All genes",]$Freq-orthocast[,2],positive=orthocast %>% select(all_of(the_i)),positive_no=tot_genes[tot_genes$cluster == the_i,]$Freq-( orthocast %>% select(all_of(the_i))))
  #another loop! for all of the ortholog types in that cluster, do a chi square to compare with freq in whole genome
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
		#the number of cluster-specific genes in that category
		#the number of cluster-specific trajectory genes not in that category
		#the fraction of genes in the genome in that category
		#the fraction of cluster-specific genes in that category
		#the chi square test statistic
		#the chisquare p value
		#the cluster of interest
    row_to_add <- data.frame(orthotype=orthoclass,all_genes_with=i_df[j,2],all_genes_without=i_df[j,3], test_genes_with=i_df[j,4],test_genes_without=i_df[j,5], fra_all_genes_with=i_df[j,2]/sum(i_df[,2]), fra_positive_genes_with=i_df[j,4]/sum(i_df[,4]), chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value,cluster=the_i)
    #add to df and repeat
    testfun <- rbind(testfun,row_to_add)
  }

}


#adjust p for multiple testing
testfun$p.adjust <- p.adjust(testfun$chisq_p,method="BH")
#write data to file

write.table(testfun,"H_clusters_ortholog_type_chi_square.tsv",quote=FALSE,row.names = FALSE,sep="\t")
	#this is supplemental table 20


#okay, now interpro domains


	#inopinata gene-domain info; generated in interproscan/workflow_notes.sh , https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/interproscan/workflow_notes.sh


inopinata_interpro <- read.table("inopinata_04.tsv", sep="\t")

#merge with expression data
names(inopinata_interpro)[names(inopinata_interpro) == 'V1'] <- 'gene'
names(inopinata_interpro)[names(inopinata_interpro) == 'V2'] <- 'interpro_id'
names(inopinata_interpro)[names(inopinata_interpro) == 'V3'] <- 'interpro_description'

big_df_merge <- merge(inopinata_interpro,ortho_merge)

#count the domains

inopinata_domain_counts <- as.data.frame(table(inopinata_interpro$interpro_id))

#get just id's and descriptions

interpro_domains_descriptions <- data.frame(interpro_id=inopinata_interpro$interpro_id,interpro_description=inopinata_interpro$interpro_description)
nrow(interpro_domains_descriptions)
#[1] 31385

interpro_domains_descriptions <- interpro_domains_descriptions[!duplicated(interpro_domains_descriptions), ]

nrow(interpro_domains_descriptions)
#[1] 6970

#ok, join descriptions and get some frequencies, differences

names(inopinata_domain_counts)[names(inopinata_domain_counts) == 'Var1'] <- 'interpro_id'
names(inopinata_domain_counts)[names(inopinata_domain_counts) == 'Freq'] <- 'count'
#merge domain counts with descriptions
inopinata_domain_counts <- merge(inopinata_domain_counts,interpro_domains_descriptions)

#get total number of gene-domain pairs
inopinata_domain_counts$total <- sum(inopinata_domain_counts$count)
#get number of gene-domain pairs NOT with a give domain
inopinata_domain_counts$num_not_domain <- inopinata_domain_counts$total-inopinata_domain_counts$count
#get fraction of gene-domain pairs among the total
inopinata_domain_counts$fra_domain <- inopinata_domain_counts$count/inopinata_domain_counts$total

#okay, domain counts for interesting clusters
#the clusters with the most interesting dynamics are 20,15,8 (increasing); 19 (max L4); and 16, 17 (decreasing)
#all genes
inopinata_domain_counts$cluster <- "All genes"
#separate by cluster
cl20 <- big_df_merge[big_df_merge$clust_k_20 == "20",]
cl15 <- big_df_merge[big_df_merge$clust_k_20 == "15",]
cl8 <- big_df_merge[big_df_merge$clust_k_20 == "8",]
cl19 <- big_df_merge[big_df_merge$clust_k_20 == "19",]
cl16 <- big_df_merge[big_df_merge$clust_k_20 == "16",]
cl17 <- big_df_merge[big_df_merge$clust_k_20 == "17",]
#count domains
cl20_ortho_cat <- as.data.frame(table(cl20$interpro_id))
cl15_ortho_cat <- as.data.frame(table(cl15$interpro_id))
cl8_ortho_cat <- as.data.frame(table(cl8$interpro_id))
cl19_ortho_cat <- as.data.frame(table(cl19$interpro_id))
cl16_ortho_cat <- as.data.frame(table(cl16$interpro_id))
cl17_ortho_cat <- as.data.frame(table(cl17$interpro_id))
#give cluster ids
cl20_ortho_cat$cluster <- "20"
cl15_ortho_cat$cluster <- "15"
cl8_ortho_cat$cluster <- "8"
cl19_ortho_cat$cluster <- "19"
cl16_ortho_cat$cluster <- "16"
cl17_ortho_cat$cluster <- "17"
#combine dfs
domain_cat_df <- rbind(cl20_ortho_cat,cl15_ortho_cat,cl8_ortho_cat,cl19_ortho_cat,cl16_ortho_cat,cl17_ortho_cat)


#finishing plot

meltdat_merge_agg_relevant_clusters <- subset(meltdat_merge_agg, cluster %in% domain_cat_df$cluster)

meltdat_merge_agg_relevant_clusters$cluster <- droplevels(meltdat_merge_agg_relevant_clusters$cluster)


meltdat_merge_agg_relevant_clusters$cluster <- factor(meltdat_merge_agg_relevant_clusters$cluster,levels = c("8","15","20","19","16","17"))


panel_a <- ggplot(meltdat_merge_agg_relevant_clusters, aes(x=stage, y=expression,group=gene)) + geom_line(alpha=0.15,size=0.15) + facet_rep_grid(~cluster) + scale_color_brewer(palette="Set1") + theme_cowplot() +ylab("Transformed expression") + xlab("Stage") + ggtitle("a")  + theme(plot.title = element_text(face = "bold",hjust = -0.015,size=24),strip.background=element_blank())


#this is supplemental figure 7A
#panel b is above but also here:


panel_b <- ggplot(ortho_cat_df, aes(x = cluster, y = Freq, fill = Var1)) + geom_col(position = "fill") + scale_fill_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="bottom",legend.title=element_blank(),legend.justification = "center",plot.title = element_text(face = "bold",hjust = -0.015,size=24)) + annotate("text",x=1, y=0.05, label=tot_genes$Freq[1],size=5) + annotate("text",x=2, y=0.05, label=tot_genes$Freq[2],size=5) + annotate("text",x=3, y=0.05, label=tot_genes$Freq[3],size=5) + annotate("text",x=4, y=0.05, label=tot_genes$Freq[4],size=5) + annotate("text",x=5, y=0.05, label=tot_genes$Freq[5],size=5) + annotate("text",x=6, y=0.05, label=tot_genes$Freq[6],size=5) + annotate("text",x=7, y=0.05, label=tot_genes$Freq[7],size=5) + xlab("Cluster") + ylab("Proportion") + ggtitle("b") + guides(fill=guide_legend(nrow=2,byrow=TRUE))

library(patchwork)
panel_a/panel_b
	#this is supplemental figure 7


#getting real tired of commenting!


#okay, domain chi square stats

#stats, domains
inopinata_domain_counts_ii <- as.data.frame(table(inopinata_interpro$interpro_id))
inopinata_domain_counts_ii$cluster <- "All_genes"

big_df_merge$interpro_id.clust <- paste(big_df_merge$interpro_id, big_df_merge$clust_k_20)
domain_cat_df <- as.data.frame(table(big_df_merge$interpro_id.clust))

domain_cat_df<- separate(domain_cat_df, Var1, into = c("Var1","cluster"), sep = " (?=[^ ]+$)")



domain_cat_df[domain_cat_df$cluster=="1",]$cluster <- "C01"
domain_cat_df[domain_cat_df$cluster=="2",]$cluster <- "C02"
domain_cat_df[domain_cat_df$cluster=="3",]$cluster <- "C03"
domain_cat_df[domain_cat_df$cluster=="4",]$cluster <- "C04"
domain_cat_df[domain_cat_df$cluster=="5",]$cluster <- "C05"
domain_cat_df[domain_cat_df$cluster=="6",]$cluster <- "C06"
domain_cat_df[domain_cat_df$cluster=="7",]$cluster <- "C07"
domain_cat_df[domain_cat_df$cluster=="8",]$cluster <- "C08"
domain_cat_df[domain_cat_df$cluster=="9",]$cluster <- "C09"
domain_cat_df[domain_cat_df$cluster=="10",]$cluster <- "C10"
domain_cat_df[domain_cat_df$cluster=="11",]$cluster <- "C11"
domain_cat_df[domain_cat_df$cluster=="12",]$cluster <- "C12"
domain_cat_df[domain_cat_df$cluster=="13",]$cluster <- "C13"
domain_cat_df[domain_cat_df$cluster=="14",]$cluster <- "C14"
domain_cat_df[domain_cat_df$cluster=="15",]$cluster <- "C15"
domain_cat_df[domain_cat_df$cluster=="16",]$cluster <- "C16"
domain_cat_df[domain_cat_df$cluster=="17",]$cluster <- "C17"
domain_cat_df[domain_cat_df$cluster=="18",]$cluster <- "C18"
domain_cat_df[domain_cat_df$cluster=="19",]$cluster <- "C19"
domain_cat_df[domain_cat_df$cluster=="20",]$cluster <- "C20"

domain_cat_df_ii <- data.frame(Var1=domain_cat_df$Var1,Freq=domain_cat_df$Freq,cluster=domain_cat_df$cluster)

domain_cat_df_iii <- rbind(inopinata_domain_counts_ii,domain_cat_df_ii)

domain_cat_df_iii$cluster <- as.factor(domain_cat_df_iii$cluster)


tot_genes <- aggregate(Freq ~ cluster, FUN=sum, data=domain_cat_df_iii)


orthocast <- dcast(domain_cat_df_iii,Var1~cluster,value.var="Freq")

orthocast[is.na(orthocast)] <- 0



testfun<-NULL

#a big loop for all of the chi square tests for domains, see above for more detailed notes on how this works.

for(i in levels(domain_cat_df_iii$cluster)[-1]){
  the_i <- i
  
  i_df <-NULL
  
  i_df <- data.frame(orthotype=orthocast$Var1, All_genes=orthocast[,2],All_genes_no= tot_genes[tot_genes$cluster == "All_genes",]$Freq-orthocast[,2],positive=orthocast %>% select(all_of(the_i)),positive_no=tot_genes[tot_genes$cluster == the_i,]$Freq-(orthocast %>% select(all_of(the_i))))

  for(j in 1:nrow(i_df)){
    orthoclass <- i_df$orthotype[j]
    the_matrix <- matrix(c(i_df[j,2],i_df[j,3],i_df[j,4],i_df[j,5]),ncol=2)
    matrixtest <- chisq.test(the_matrix)
    row_to_add <- data.frame(orthotype=orthoclass,all_genes_with=i_df[j,2],all_genes_without=i_df[j,3], test_genes_with=i_df[j,4],test_genes_without=i_df[j,5], fra_all_genes_with=i_df[j,2]/sum(i_df[,2]), fra_positive_genes_with=i_df[j,4]/sum(i_df[,4]), chisq_statistic=matrixtest$statistic,chisq_df=matrixtest$parameter,chisq_p=matrixtest$p.value,cluster=the_i)
    testfun <- rbind(testfun,row_to_add)
  }

}
	#this may take a while

#adjust for p values
testfun$p.adjust <- p.adjust(testfun$chisq_p,method="BH")
#replace the id
names(testfun)[names(testfun) == 'orthotype'] <- 'interpro_id'
#not in the paper, but on github
write.table(testfun,"H_clusters_domains_chi_square.tsv",quote=FALSE,row.names = FALSE,sep="\t")

#just get significant domains

sig_ids <- testfun[testfun$p.adjust < 0.05,]

#merge with descriptions

sig_ids_merge <- merge(sig_ids,interpro_domains_descriptions)
#write to a file
write.table(sig_ids_merge,"hclust_sig_interpro_ids_with_descriptions.tsv",quote=FALSE,row.names = FALSE,sep="\t")
	#this is supplemental table 19
