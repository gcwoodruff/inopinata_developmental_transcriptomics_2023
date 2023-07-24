library(airway)
library(tximport)
library(DESeq2)
library(PoiClaClu)
library(reshape)

setwd("/Users/gavin/genome/inopinata_RNAseq_1-22-20/MS_code_7-2023/")

  #put directory with the "quants" folder here
dir <- "/Users/gavin/genome/inopinata_RNAseq_1-22-20/MS_code_7-2023/"

	#file

coldata <- read.table("sample_table_b.tsv", sep="\t", header=T)

coldata$files <- file.path(dir, "quants", coldata$sample_id)

files <- file.path(dir, "quants", coldata$sample_id)

tx2gene <- read.table("tx2gene.tsv", sep="\t", header=T)

txi <- tximport(coldata$files, type = "salmon", tx2gene = tx2gene)



txi$group <- factor(paste0(txi$Species, txi$Stage))
coldata$group <- factor(paste0(coldata$Species, coldata$Stage))

dat <- DESeqDataSetFromTximport(txi, coldata, design= ~ group)

countdata <- round(assays(dat)[["counts"]])

nrow(dat)
#[1] 10718

keep <- rowSums(counts(dat)) > 1
dat <- dat[keep,]
nrow(dat)
#[1] 10717 ; we lost one!


#data transformation

rld <- rlog(dat, blind = FALSE)
head(assay(rld,assaywithDimnames=TRUE), 3)

#              1        2        3        4        5        6        7        8
#2L52.1 7.708664 7.870636 7.997748 7.854591 7.813753 6.948719 6.791148 7.263215
#4R79.2 4.388329 3.449245 3.190055 3.997312 3.758136 4.829324 4.950777 5.070108
#aagr-1 8.454997 7.127890 8.285564 7.577651 7.774986 8.469207 8.869573 8.779053
#              9       10       11       12       13       14       15       16
#2L52.1 7.093394 7.331547 7.604789 6.958331 7.210628 7.316678 6.451683 6.519520
#4R79.2 5.033134 5.034328 4.026636 4.903430 3.915434 4.058440 4.804630 5.223304
#aagr-1 8.680929 8.927315 8.732071 8.076126 8.455288 8.569642 9.015869 9.053834
#             17       18        19        20       21        22        23
#2L52.1 6.646279 6.332676  6.274892  6.677204 6.392063  6.077938  6.531338
#4R79.2 5.479098 4.499062  3.561607  4.129477 5.210211  4.936470  5.069274
#aagr-1 8.970998 8.618756 10.526390 10.366639 9.976514 10.039445 10.056347
#              24        25        26       27        28
#2L52.1  6.729144  6.855788  6.938440 7.176714  6.909910
#4R79.2  4.592185  4.238838  4.315236 3.798145  3.693291
#aagr-1 10.136719 10.141410 10.171365 9.801042 10.191233




#sample distances

sampleDists <- dist(t(assay(rld)))
sampleDists




poisd <- PoissonDistance(t(counts(dat)))


samplePoisDistMatrix <- as.matrix( poisd$dd )


rownames(samplePoisDistMatrix) <- dat@colData$sample_id
colnames(samplePoisDistMatrix) <- dat@colData$sample_id

#write these distances to a file
#write.table(as.data.frame(samplePoisDistMatrix), file="sample_pois_dist_matrix.tsv", sep="\t", quote = FALSE)

m <- as.matrix(sampleDists)
m
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2) <- c("c1", "c2", "distance")
m2

#write these distances to a file
#write.table(as.data.frame(m2), file= 'sample_pairwise_distances.tsv',sep="\t", row.names = FALSE, quote = FALSE)


#differentatial expression between species at each stage


DS2_dat <- DESeq(dat)

resultsNames(DS2_dat)

res_L3 <- results(DS2_dat, contrast=c("group","inopinataL3","elegansL3"))
res_L4 <- results(DS2_dat, contrast=c("group","inopinataL4","elegansL4"))
res_adult <- results(DS2_dat, contrast=c("group","inopinataadult","elegansadult"))


#what are the number of significant differentially expressed genes between species at each stage?
sum(res_L3$padj < 0.05, na.rm=TRUE)
#[1] 6057
sum(res_L4$padj < 0.05, na.rm=TRUE)
#[1] 5900
sum(res_adult$padj < 0.05, na.rm=TRUE)
#[1] 7126

#the fraction?
sum(res_L3$padj < 0.05, na.rm=TRUE)/10717
#[1] 0.5651768

sum(res_L4$padj < 0.05, na.rm=TRUE)/10717
#[1] 0.5505272

sum(res_adult$padj < 0.05, na.rm=TRUE)/10717
#[1] 0.6649249

#print the differential expression analysis for each stage to a tsv

	#these are results used for analysis and figures!
write.table(as.data.frame(res_L3), file= 'L3_stats_results.tsv',sep="\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(res_L4), file= 'L4_stats_results.tsv',sep="\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(res_adult), file= 'adult_stats_results.tsv',sep="\t", row.names = FALSE, quote = FALSE)

#print all the regularized log transformed data to a file
	#these are results used for analysis and figures!
write.table(assay(rld), file="rlog_transformed_transcript_counts.tsv", sep="\t", quote = FALSE)




#now include an interaction in the model

#include the interaction
ddsTC <- DESeqDataSet(dat, ~ Species + Stage + Species:Stage)

#do the LRT with the reduced model
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ Species + Stage)
#get the results
resTC <- results(ddsTC)
#this gets the Species:Adult-L4 interaction
resTC_iL4 <- results(ddsTC,name=c("Speciesinopinata.StageL4"))

#check the data (here, although DeSeq2 calls the column "log2FoldChange", here it is the interaction term)
head(resTC_iL4[order(-resTC_iL4$log2FoldChange),],)

#confirm this-- that column has the interaction coefficients
ddsTC_coef <- coef(ddsTC)

head(resTC_iL4)
head(ddsTC_coef)

	#these _are_ the interaction coefficients, great!

#write the data to a file ; this is data used for analysis and figures
write.table(resTC_iL4[order(-resTC_iL4$log2FoldChange),], file= 'time_course_LRT_results_Speciesinopinata.StageL4_by_int_term.tsv',sep="\t", quote = FALSE)





#order by asbolute value of interaction term 
dat_order_abs_val_int_h_l <- dat_sig[order(-dat_sig$int_term_abs),]

dat_order_abs_val_int_h_l_list <- head(dat_order_abs_val_int_h_l,720)

write.table(dat_order_abs_val_int_h_l_list, file= 'LIST_A_top_ten_perc_abs_val_Species--L4-adult_int_term_only_sig.tsv',sep="\t", quote = FALSE)


#order by value of interaction term (high to low; top positive interactions) ; this is the positive interactions list


dat_order_int_h_l <- dat_sig[order(-dat_sig$int_term_x_neg_one),]

dat_order_int_h_l_list <- head(dat_order_int_h_l,720)

write.table(dat_order_int_h_l_list, file= 'LIST_B_top_ten_perc_positive_Species--L4-adult_int_term_only_sig.tsv',sep="\t", quote = FALSE)


#order by value of interaction term (low to high; top negative interactions)  ; this is the positive interactions list


dat_order_int_l_h <- dat_sig[order(dat_sig$int_term_x_neg_one),]

dat_order_int_l_h_list <- head(dat_order_int_l_h,720)

write.table(dat_order_int_l_h_list, file= 'LIST_C_top_ten_perc_negative_Species--L4-adult_int_term_only_sig.tsv',sep="\t", quote = FALSE)



#order by p value


dat_order_int_p <- dat_sig[order(-dat_sig$padj),]

dat_order_int_p_list <- head(dat_order_int_p,720)

write.table(dat_order_int_p_list, file= 'LIST_D_top_ten_perc_p_value_Species--L4-adult_int_term_only_sig.tsv',sep="\t", quote = FALSE)


#PCA



pca <- prcomp(t(assay(rld)))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

intgroup.df <- as.data.frame(colData(rld)[, c( "Species", "Stage"), drop=FALSE])
  # add the intgroup factors together to create a new grouping factor
group <- if (length(c( "Species", "Stage")) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(rld)[[c( "Species", "Stage")]]
}
#PC1 and PC2 percent variance explained
round(percentVar[1] * 100)
round(percentVar[2] * 100)

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(rld))


	#used for analysis and plotting
write.table(d, file= 'all_genes_pca.tsv',sep="\t", row.names = FALSE, quote = FALSE)

