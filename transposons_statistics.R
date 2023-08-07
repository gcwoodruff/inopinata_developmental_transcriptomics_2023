library(reshape2)
setwd("/Users/gavin/genome/inopinata_RNAseq_1-22-20/figures_MS_7-2023")


trans_dat <- read.table("rlog_transformed_transcript_counts_transposon_tags.tsv", header=TRUE,sep="\t")



trans_dat <- melt(trans_dat)

names(trans_dat)[names(trans_dat) == 'variable'] <- 'sample_id'
names(trans_dat)[names(trans_dat) == 'value'] <- 'transcript_count'

trans_dat$stage[trans_dat$sample_id == "sample_1"] <- "L3"
trans_dat$stage <- as.factor(trans_dat$stage)
levels(trans_dat$stage) <- c("L3","L4","adult")

trans_dat$stage[trans_dat$sample_id == "sample_1"] <- "L3"
trans_dat$stage[trans_dat$sample_id == "sample_2"] <- "L3"
trans_dat$stage[trans_dat$sample_id == "sample_3"] <- "L3"
trans_dat$stage[trans_dat$sample_id == "sample_4"] <- "L3"
trans_dat$stage[trans_dat$sample_id == "sample_5"] <- "L3"
trans_dat$stage[trans_dat$sample_id == "sample_6"] <- "L4"
trans_dat$stage[trans_dat$sample_id == "sample_7"] <- "L4"
trans_dat$stage[trans_dat$sample_id == "sample_9"] <- "L4"
trans_dat$stage[trans_dat$sample_id == "sample_10"] <- "L4"
trans_dat$stage[trans_dat$sample_id == "sample_11"] <- "adult"
trans_dat$stage[trans_dat$sample_id == "sample_12"] <- "adult"
trans_dat$stage[trans_dat$sample_id == "sample_13"] <- "adult"
trans_dat$stage[trans_dat$sample_id == "sample_14"] <- "adult"
trans_dat$stage[trans_dat$sample_id == "sample_15"] <- "adult"

trans_dat$transposon_status <- as.factor(trans_dat$transposon_status)
levels(trans_dat$transposon_status)[levels(trans_dat$transposon_status)=="NOT transposon-aligning"] <- "No"
levels(trans_dat$transposon_status)[levels(trans_dat$transposon_status)=="transposon-aligning"] <- "Yes"


trans_dat$gene.stage.transposon <- as.factor(paste(trans_dat$gene_id, trans_dat$stage, trans_dat$transposon_status))


new_trans_dat_df <- aggregate(transcript_count ~ gene.stage.transposon, FUN=mean, data=trans_dat)



library(tidyr)
new_trans_dat_df <- separate(new_trans_dat_df, gene.stage.transposon, into = c("gene","stage","transposon"), sep = " ")


new_trans_dat_df$transposon <- as.factor(new_trans_dat_df$transposon)

new_trans_dat_df$stage <- as.factor(new_trans_dat_df$stage)


new_trans_dat_df$transposon <- factor(new_trans_dat_df$transposon, levels=c("Yes","No"))

new_trans_dat_df$stage <- factor(new_trans_dat_df$stage, levels=c("L3", "L4", "adult"))
levels(new_trans_dat_df$stage)[levels(new_trans_dat_df$stage)=="adult"] <- "Adult"


library(effsize)
L3tx <- new_trans_dat_df[new_trans_dat_df$stage == "L3",]
L4tx <- new_trans_dat_df[new_trans_dat_df$stage == "L4",]
Adulttx <- new_trans_dat_df[new_trans_dat_df$stage == "Adult",]

perc_diff = function(w,z){
	((w-z)/z)*100
}

perc_diff2 = function(w,z){
	(1-(z-w)/z)*100
}

#L3
#summary
summary(L3tx[L3tx$transposon == "No",]$transcript_count)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -2.042   1.320   5.152   4.892   8.236  17.309

summary(L3tx[L3tx$transposon == "Yes",]$transcript_count)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-2.0425 -0.6973  1.0323  1.9724  3.9173 15.3671

#% difference
perc_diff(mean(L3tx[L3tx$transposon == "Yes",]$transcript_count),mean(L3tx[L3tx$transposon == "No",]$transcript_count))
#[1] -59.68437

#cohen's d effect size
cohen.d(L3tx[L3tx$transposon == "Yes",]$transcript_count,L3tx[L3tx$transposon == "No",]$transcript_count)
#Cohen's d
#
#d estimate: -0.7291488 (medium)
#95 percent confidence interval:
#     lower      upper
#-0.7994007 -0.6588968
#

##wilcoxon rank sum test

wilcox.test(L3tx[L3tx$transposon == "Yes",]$transcript_count,L3tx[L3tx$transposon == "No",]$transcript_count,exact = FALSE)
#	Wilcoxon rank sum test with continuity correction
#
#data:  L3tx[L3tx$transposon == "Yes", ]$transcript_count and L3tx[L3tx$transposon == "No", ]$transcript_count
#W = 4477356, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0





#L4
#summary
summary(L4tx[L4tx$transposon == "No",]$transcript_count)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -2.173   2.569   6.192   5.467   8.442  16.557

summary(L4tx[L4tx$transposon == "Yes",]$transcript_count)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-2.1731 -0.6901  1.3876  2.3485  4.7711 14.6017

#% difference
perc_diff(mean(L4tx[L4tx$transposon == "Yes",]$transcript_count),mean(L4tx[L4tx$transposon == "No",]$transcript_count))
#[1] -57.04592

#cohen's d effect size
cohen.d(L4tx[L4tx$transposon == "Yes",]$transcript_count,L4tx[L4tx$transposon == "No",]$transcript_count)
#Cohen's d
#
#d estimate: -0.8100293 (large)
#95 percent confidence interval:
#     lower      upper
#-0.8803677 -0.7396908

##wilcoxon rank sum test

wilcox.test(L4tx[L4tx$transposon == "Yes",]$transcript_count,L4tx[L4tx$transposon == "No",]$transcript_count,exact = FALSE)
#	Wilcoxon rank sum test with continuity correction
#
#data:  L4tx[L4tx$transposon == "Yes", ]$transcript_count and L4tx[L4tx$transposon == "No", ]$transcript_count
#W = 4298340, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0








#Adult
#summary
summary(Adulttx[Adulttx$transposon == "No",]$transcript_count)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -2.093   2.309   6.003   5.401   8.621  18.117

summary(Adulttx[Adulttx$transposon == "Yes",]$transcript_count)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-2.0925 -0.7459  1.1009  2.2112  4.5039 14.0341

#% difference
perc_diff(mean(Adulttx[Adulttx$transposon == "Yes",]$transcript_count),mean(Adulttx[Adulttx$transposon == "No",]$transcript_count))
#[1] -59.06062

#cohen's d effect size
cohen.d(Adulttx[Adulttx$transposon == "Yes",]$transcript_count,Adulttx[Adulttx$transposon == "No",]$transcript_count)
#Cohen's d
#
#d estimate: -0.8177325 (large)
#95 percent confidence interval:
#     lower      upper
#-0.8880796 -0.7473853

##wilcoxon rank sum test

wilcox.test(Adulttx[Adulttx$transposon == "Yes",]$transcript_count,Adulttx[Adulttx$transposon == "No",]$transcript_count,exact = FALSE)
#	Wilcoxon rank sum test with continuity correction
#
#data:  Adulttx[Adulttx$transposon == "Yes", ]$transcript_count and Adulttx[Adulttx$transposon == "No", ]$transcript_count
#W = 4269414, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0










