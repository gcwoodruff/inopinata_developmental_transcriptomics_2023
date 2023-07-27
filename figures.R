
library(ggplot2)
library(ggforce)
library(cowplot)
library(patchwork)
library(reshape2)
library(lemon)
library(GGally)

setwd("/Users/gavin/genome/inopinata_RNAseq_1-22-20/figures_MS_7-2023")

#Figure 1
dat <- read.table("body_size.tsv", header=TRUE,sep="\t")

dat$Species <- dat$species

a <- ggplot(dat, aes(x = day, y = length)) + geom_rect(aes(xmin=2.75,xmax=4.25,ymin=0,ymax=2000),fill="khaki1",alpha=0.01) + geom_jitter(aes(colour=Species),size=1,width = 0.12,alpha=0.5)  + stat_summary(aes(group=Species,colour=Species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25,size=0.25) + stat_summary(aes(group=Species,colour=Species,linetype=Species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "line", width = 0.25) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.text = element_text(face = "italic")) + xlab("Days since embryo laid") + ylab("Length (microns)") + scale_y_continuous(limits=c(0,2000),breaks=c(0,500,1000,1500,2000)) + scale_x_continuous(limits=c(0.5,7.5),breaks=c(1:7))


stagedat <- read.table("length_width_data_2022_revision.tsv", header=TRUE,sep="\t",stringsAsFactors = TRUE)

gcwdat <- stagedat[stagedat$observer_year == "GCW Evol. Lett. 2018",]

L3 <- gcwdat[gcwdat$stage == "L3",]
L4 <- gcwdat[gcwdat$stage == "L4",]
adult <- gcwdat[gcwdat$stage == "Gravid adult",]

plotdat <- rbind(L3,L4,adult)

plotdat$stage <- droplevels(plotdat$stage)

plotdat$stage <- factor(plotdat$stage, levels=c("L3", "L4", "Gravid adult"))
levels(plotdat$stage)[levels(plotdat$stage)=="Gravid adult"] <- "Adult"

plotdat$species <- droplevels(plotdat$species)
plotdat$Species <- plotdat$species
plotdat$Stage <- plotdat$stage

b <- ggplot(plotdat, aes(x = Stage, y = length)) + geom_sina(aes(colour=Species),size=1,scale="width") + stat_summary(aes(group=Species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1",name="Species") + theme_cowplot() + theme(legend.text = element_text(face = "italic")) + ylab("Length (microns)") + scale_y_continuous(limits=c(0,2000),breaks=c(0,500,1000,1500,2000)) 

#this is figure 1
a+b

#ggsave("figure_1_body_size.png",bg="white", units="in", width= 6,height=4.5)
ggsave("figure_1_body_size.pdf",useDingbats=FALSE, units="in",  width= 6,height=4.5)


a <- ggplot(dat, aes(x = day, y = length)) + geom_rect(aes(xmin=2.75,xmax=4.25,ymin=0,ymax=2000),fill="khaki1",alpha=0.01) + geom_jitter(aes(colour=Species),size=1,width = 0.12)  + stat_summary(aes(group=Species,colour=Species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25,size=0.25) + stat_summary(aes(group=Species,colour=Species,linetype=Species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "line", width = 0.25) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.position="none") + xlab("Days since embryo laid") + ylab("Length (microns)") + scale_y_continuous(limits=c(0,2000),breaks=c(0,500,1000,1500,2000)) + scale_x_continuous(limits=c(0.5,7.5),breaks=c(1:7))
b <- ggplot(plotdat, aes(x = Stage, y = length)) + geom_sina(aes(colour=Species),size=1,scale="width") + stat_summary(aes(group=Species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1",name="Species") + theme_cowplot() + theme(legend.position="none") + ylab("Length (microns)") + scale_y_continuous(limits=c(0,2000),breaks=c(0,500,1000,1500,2000))

#figure 1
a+b

ggsave("figure_1_body_size_no_legend.pdf",useDingbats=FALSE, units="in",  width= 6,height=4.5)



#figure 2


#volcano plots smbe poster

L3dat <- read.table("L3_stats_results.tsv", sep="\t", header=TRUE)

L3dat$neg.log.p.adj <- -(log10(L3dat$padj))

L3dat$is.sig <- L3dat$padj < 0.05

L3dat.om <- na.omit(L3dat)

L3dat.om$is.sig <- as.factor(L3dat.om$is.sig)

levels(L3dat.om$is.sig)[levels(L3dat.om$is.sig)=="FALSE"] <- "No"
levels(L3dat.om$is.sig)[levels(L3dat.om$is.sig)=="TRUE"] <- "Yes"

a <-ggplot(L3dat.om, aes(x = log2FoldChange, y = neg.log.p.adj)) + geom_point(size=0.3,alpha=0.5, aes(colour=is.sig)) + scale_colour_manual(name= "p < 0.05?",values=c("black","red")) + theme_cowplot() + theme(legend.position="none",plot.title = element_text(face = "plain")) + xlab(expression( ~Log[2]~paste("fold change"))) + ylab(expression( ~-Log[10]~paste(italic(" P")))) + ggtitle("L3") + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))



L4dat <- read.table("L4_stats_results.tsv", sep="\t", header=TRUE)

L4dat$neg.log.p.adj <- -(log10(L4dat$padj))

L4dat$is.sig <- L4dat$padj < 0.05

L4dat.om <- na.omit(L4dat)


L4dat.om$is.sig <- as.factor(L4dat.om$is.sig)

levels(L4dat.om$is.sig)[levels(L4dat.om$is.sig)=="FALSE"] <- "No"
levels(L4dat.om$is.sig)[levels(L4dat.om$is.sig)=="TRUE"] <- "Yes"


b <- ggplot(L4dat.om, aes(x = log2FoldChange, y = neg.log.p.adj)) + geom_point(size=0.3,alpha=0.5, aes(colour=is.sig)) + scale_colour_manual(name= "p < 0.05?",values=c("black","red")) + theme_cowplot()  + xlab(expression( ~Log[2]~paste("fold change"))) +theme(axis.title.y=element_blank(),plot.title = element_text(face = "plain")) + ylab(expression( ~-Log[10]~paste(italic(" P")))) + ggtitle("L4") + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))



adultdat <- read.table("adult_stats_results.tsv", sep="\t", header=TRUE)

adultdat$neg.log.p.adj <- -(log10(adultdat$padj))

adultdat$is.sig <- adultdat$padj < 0.05

adultdat.om <- na.omit(adultdat)

adultdat.om$is.sig <- as.factor(adultdat.om$is.sig)

levels(adultdat.om$is.sig)[levels(adultdat.om$is.sig)=="FALSE"] <- "No"
levels(adultdat.om$is.sig)[levels(adultdat.om$is.sig)=="TRUE"] <- "Yes"


c <- ggplot(adultdat.om, aes(x = log2FoldChange, y = neg.log.p.adj)) + geom_point(size=0.3,alpha=0.5, aes(colour=is.sig)) + scale_colour_manual(name= "p < 0.05?",values=c("black","red"))  + theme_cowplot() +theme(legend.position="none",plot.title = element_text(face = "plain")) + xlab(expression( ~Log[2]~paste("fold change"))) + ylab(expression( ~-Log[10]~paste(italic(" P")))) + ggtitle("Adult") + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))



intdat <- read.table("time_course_LRT_results_Speciesinopinata.StageL4_by_int_term.tsv", sep="\t", header=TRUE)
intdat$int_term_x_neg_one <- intdat$log2FoldChange*(-1)
intdat$neg.log.p.adj <- -(log10(intdat$padj))
intdat$int_term_abs <- abs(intdat$log2FoldChange)
intdat$is.sig <- intdat$padj < 0.05
intdat.om <- na.omit(intdat)
intdat.om$is.sig <- as.factor(intdat.om$is.sig)
levels(intdat.om$is.sig)[levels(intdat.om$is.sig)=="FALSE"] <- "No"
levels(intdat.om$is.sig)[levels(intdat.om$is.sig)=="TRUE"] <- "Yes"


d <- ggplot(intdat.om, aes(x = int_term_x_neg_one, y = neg.log.p.adj)) + geom_point(size=0.3,alpha=0.25, aes(colour=is.sig)) + scale_colour_manual(values=c("black","red")) + theme_cowplot() +theme(axis.title.y=element_blank(),legend.position="none",plot.title = element_text(face = "plain")) + xlab("Species:L4-adult interaction coefficient") + ylab(bquote(~-log[10] ~ italic(p)))+ labs(colour=bquote(italic(p)<0.05)) + ggtitle("L4:Adult Interaction")

#this is figure 2
(a+b)/(c+d)


ggsave("figure_2.pdf",useDingbats=FALSE, units="in",  width= 7.5,height=5.5)


#figure 3


dat <- read.table("time_course_LRT_results_Speciesinopinata.StageL4_by_int_term.tsv", sep="\t", header=TRUE)

dat$int_term_x_neg_one <- dat$log2FoldChange*(-1)
dat$neg.log.p.adj <- -(log10(dat$padj))
dat$int_term_abs <- abs(dat$log2FoldChange)
dat$is.sig <- dat$padj < 0.05


dat_sig <- subset(dat, padj < 0.05)

dat_order_int_pos <- dat_sig[order(-dat_sig$int_term_x_neg_one),]
dat_order_int_pos_top_ten <- head(dat_order_int_pos,10)


dat_order_int_neg <- dat_sig[order(dat_sig$int_term_x_neg_one),]
dat_order_int_neg_top_ten <- head(dat_order_int_neg,10)



rlogdat <- read.table("rlog_transformed_transcript_counts.tsv", sep="\t", header=TRUE)

rlogdat <- melt(rlogdat)

names(rlogdat)[names(rlogdat) == 'variable'] <- 'sample_id'
names(rlogdat)[names(rlogdat) == 'value'] <- 'transcript_count'


rlogdat$species[rlogdat$sample_id == "s1"] <- "C. inopinata"
rlogdat$species <- as.factor(rlogdat$species)
levels(rlogdat$species) <- c("C. elegans", "C. inopinata")


rlogdat$species[rlogdat$sample_id == "s1"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s2"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s3"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s4"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s5"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s6"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s7"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s9"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s10"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s11"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s12"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s13"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s14"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s15"] <- "C. inopinata"
rlogdat$species[rlogdat$sample_id == "s16"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s17"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s18"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s19"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s20"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s21"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s22"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s23"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s24"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s25"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s26"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s27"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s28"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s29"] <- "C. elegans"
rlogdat$species[rlogdat$sample_id == "s30"] <- "C. elegans"




rlogdat$stage[rlogdat$sample_id == "s1"] <- "L3"
rlogdat$stage <- as.factor(rlogdat$stage)
levels(rlogdat$stage) <- c("L3","L4","adult")




rlogdat$stage[rlogdat$sample_id == "s1"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s2"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s3"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s4"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s5"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s6"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s7"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s9"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s10"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s11"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s12"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s13"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s14"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s15"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s16"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s17"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s18"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s19"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s20"] <- "L3"
rlogdat$stage[rlogdat$sample_id == "s21"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s22"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s23"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s24"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s25"] <- "L4"
rlogdat$stage[rlogdat$sample_id == "s26"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s27"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s28"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s29"] <- "adult"
rlogdat$stage[rlogdat$sample_id == "s30"] <- "adult"

rlogdat$gen_id <- as.factor(rlogdat$gen_id)

hch_1_rlogdat <-rlogdat[rlogdat$gen_id == "hch-1",]
hch_1_rlogdat <- hch_1_rlogdat[hch_1_rlogdat$stage != "L3",]
hch_1_rlogdat$gen_id <- droplevels(hch_1_rlogdat$gen_id)
hch_1_rlogdat$stage <- droplevels(hch_1_rlogdat$stage)

col_81_rlogdat <-rlogdat[rlogdat$gen_id == "col-81",]
col_81_rlogdat <- col_81_rlogdat[col_81_rlogdat$stage != "L3",]
col_81_rlogdat$gen_id <- droplevels(col_81_rlogdat$gen_id)
col_81_rlogdat$stage <- droplevels(col_81_rlogdat$stage)


a <- ggplot(dat_order_int_pos_top_ten, aes(y=reorder(gene, int_term_abs), x=int_term_abs)) + geom_bar(stat='identity',fill="#9ecae1") + theme_cowplot() +xlab("|Interaction Term|") + ylab("Gene") + theme(axis.text.x = element_text(colour="black", size=18), axis.text.y = element_text(colour="black", size=18),axis.title=element_text(size=20),title=element_text(size=20)) + ggtitle("Positive interactions")



b <- ggplot(dat_order_int_neg_top_ten, aes(y=reorder(gene, int_term_abs), x=int_term_abs)) + geom_bar(stat='identity',fill="#9ecae1") + theme_cowplot() +xlab("|Interaction Term|") + ylab("Gene") + theme(axis.text.x = element_text(colour="black", size=18), axis.text.y = element_text(colour="black", size=18),axis.title=element_text(size=20),title=element_text(size=20)) + ggtitle("Negative interactions")



c <- ggplot(hch_1_rlogdat, aes(x = stage, y = transcript_count)) + geom_sina(aes(colour=species),size=2,scale="width",position=position_dodge(width=0.75)) + stat_summary(aes(group=species,colour=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25,position = position_dodge(width = 0.75)) + stat_summary(aes(group=species,colour=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "line",position = position_dodge(width = 0.75)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.text=element_text(colour="black", size=13,face="italic"),plot.title = element_text(face="italic")) + xlab("Stage") + ylab("Transformed transcript count") + ggtitle(unique(levels(hch_1_rlogdat$gen_id)))

d <- ggplot(col_81_rlogdat, aes(x = stage, y = transcript_count)) + geom_sina(aes(colour=species),size=2,scale="width",position=position_dodge(width=0.75)) + stat_summary(aes(group=species,colour=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25,position = position_dodge(width = 0.75)) + stat_summary(aes(group=species,colour=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "line",position = position_dodge(width = 0.75)) + scale_colour_brewer(palette="Set1") + theme_cowplot() + theme(legend.text=element_text(colour="black", size=13,face="italic"),plot.title = element_text(face="italic")) + xlab("Stage") + ylab("Transformed transcript count") + ggtitle(unique(levels(col_81_rlogdat$gen_id)))

#this is figure 3
(a+c)/(b+d)


ggsave("figure_3.pdf",useDingbats=FALSE, units="in",  width= 9,height=8)



#figure 4, WB ontology


anat_dat <- read.table("anatomy_plot.csv", sep="\t", header=TRUE)
	#see get_gene_lists.R for the gene lists used to generate these data (List A, B, and C respectively)
	#Here, Lists B (positive interactions) and C (negative interactions) were used
	#These genes were used in WormBases's gene set enrichment analysis tool: https://wormbase.org/tools/enrichment/tea/tea.cgi
	#the data for gene enrichment analyses are in folder "WB_gene_enrichment_tool"

anat_dat_pos <- anat_dat[anat_dat$pos_or_neg == "Positive interaction",]
anat_dat_neg <- anat_dat[anat_dat$pos_or_neg == "Negative interaction",]



a <- ggplot(anat_dat_pos, aes(y=reorder(Term, Enrichment.Fold.Change), x=Enrichment.Fold.Change)) + geom_bar(stat='identity',fill="#9239F6") + theme_cowplot(font_size = 20) +xlab("Enrichment fold change") + ylab("Anatomy term") + ggtitle("Positive interactions")



b <- ggplot(anat_dat_neg, aes(y=reorder(Term, Enrichment.Fold.Change), x=Enrichment.Fold.Change)) + geom_bar(stat='identity',fill="#FF0076") + theme_cowplot(font_size = 20) +xlab("Enrichment fold change") + ylab("Anatomy term") + ggtitle("Negative interactions")


phen_dat <- read.table("phenotype_plot.tsv", sep="\t", header=TRUE)
	#see get_gene_lists.R for the gene lists used to generate these data (List A, B, and C respectively)
	#Here, Lists B (positive interactions) and C (negative interactions) were used
	#These genes were used in WormBases's gene set enrichment analysis tool: https://wormbase.org/tools/enrichment/tea/tea.cgi
	#the data for gene enrichment analyses are in folder "WB_gene_enrichment_tool"

phen_dat_pos <- phen_dat[phen_dat$pos_or_neg == "Positive Interaction",]
phen_dat_neg <- phen_dat[phen_dat$pos_or_neg == "Negative Interaction",]



c <- ggplot(phen_dat_pos, aes(y=reorder(Term, Enrichment.Fold.Change), x=Enrichment.Fold.Change)) + geom_bar(stat='identity',fill="#05ffa1") + theme_cowplot(font_size = 20) +xlab("Enrichment fold change") + ylab("Phenotype term") + ggtitle("Positive interactions")



d <- ggplot(phen_dat_neg, aes(y=reorder(Term, Enrichment.Fold.Change), x=Enrichment.Fold.Change)) + geom_bar(stat='identity',fill="#e5e187") + theme_cowplot(font_size = 20) +xlab("Enrichment fold change") + ylab("Phenotype term") + ggtitle("Negative interactions")


#this is figure 4
(a+c)/(b+d)



ggsave("figure_4.pdf",useDingbats=FALSE, units="in",  width= 9,height=8)


#figure 5 wormexp results


wormexp_dat <- read.table("wormexp_plot.csv", sep="\t", header=TRUE)
	#see get_gene_lists.R for the gene lists used to generate these data (List A, B, and C respectively)
	#Here, Lists B (positive interactions) and C (negative interactions) were used
	#These genes were used in the WormExp tool: https://wormexp.zoologie.uni-kiel.de/wormexp/
	#the output data from this tool are in folder "wormexp"


wormexp_dat$perc_overlap <- (wormexp_dat$Counts/wormexp_dat$ListSize)*100

wormexp_dat_pos <- wormexp_dat[wormexp_dat$pos_or_neg == "Positive interaction",]
wormexp_dat_neg <- wormexp_dat[wormexp_dat$pos_or_neg == "Negative interaction",]



a <- ggplot(wormexp_dat_pos, aes(y=reorder(Term, perc_overlap), x=perc_overlap)) + geom_bar(stat='identity',fill="#D96237") + theme_cowplot(font_size = 20) +xlab("% overlap") + ylab("Previous experiment") + ggtitle("Positive interactions")




b <- ggplot(wormexp_dat_neg, aes(y=reorder(Term, perc_overlap), x=perc_overlap)) + geom_bar(stat='identity',fill="#55FFFF") + theme_cowplot(font_size = 20) +xlab("% overlap") + ylab("Previous experiment") + ggtitle("Negative interactions")

#this is figure 5
a/b


ggsave("figure_5.pdf",units="in",width= 9,height=8,useDingbats=FALSE)




#figure 6
	#these data generated from L3_stats_results.tsv, L4_stats_results.tsv, adult_stats_results.tsv, and gumienny_savage-dunn_tgf-beta_list.txt. see get_tgf-beta_genes.sh
tgfbeta_dat <- read.table("all_stages_LRT_results_TGF-beta.tsv", sep="\t", header=TRUE)

tgfbeta_dat$stage <- factor(tgfbeta_dat$stage, levels = c("L3","L4","Adult"))



tgfbeta_dat$is.sig <- tgfbeta_dat$padj < 0.05

tgfbeta_dat$is.sig <- as.factor(tgfbeta_dat$is.sig)

levels(tgfbeta_dat$is.sig)[levels(tgfbeta_dat$is.sig)=="FALSE"] <- "No"
levels(tgfbeta_dat$is.sig)[levels(tgfbeta_dat$is.sig)=="TRUE"] <- "Yes"


tgfbeta_dat$gene_id <- factor(tgfbeta_dat$gene_id, levels = c("dbl-1","daf-7","unc-129","tig-2","tig-3","sma-6","daf-1","daf-4","sma-2","sma-3","daf-8","daf-14","sma-4","daf-3","tag-68","sma-9","lin-31","mab-31","daf-5","daf-12","lon-1","lon-2","sma-10","crm-1","drag-1","adt-2"))


cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



#this is figure 6
ggplot(tgfbeta_dat, aes(x = gene_id, y = log2FoldChange))  + geom_stripped_cols(odd = "white", even = "gray90")  + geom_col(aes(fill=stage,colour=is.sig),position=position_dodge2()) + geom_linerange(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),position=position_dodge2(0.9)) + scale_fill_manual(values=cbbPalette)  + scale_colour_manual(name= "p < 0.05?",values=c("white","black")) + theme(panel.border=element_blank(),panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour="black"), axis.text.x = element_text(colour="black", size=17,face="italic", angle = 45, hjust = 1), axis.text.y = element_text(colour="black", size=17),legend.text=element_text(colour="black", size=17),legend.title=element_text(colour="black", size=17),axis.title=element_text(size=19), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 17, colour = "black"),legend.key=element_blank())  + xlab("Gene") + ylab(expression( ~Log[2]~paste("fold change"))) + geom_segment(aes(x = 0.6, y = 6.4, xend = 5.6, yend = 6.4)) + annotate(geom="text", size = 5,x=3, y= 6.6, label="Ligands") + geom_segment(aes(x = 5.6, y = 2, xend = 8.4, yend = 2)) + annotate(geom="text", size = 5,x=7, y= 2.2, label="Receptors") + geom_segment(aes(x = 8.6, y = 3.5, xend = 15.4, yend = 3.5)) + annotate(geom="text", size = 5,x=12, y= 3.7, label="SMADS") + geom_segment(aes(x = 15.6, y = 4, xend = 20.5, yend = 4)) + annotate(geom="text", size = 5,x=18, y= 4.2, label="Transcription factors") + geom_segment(aes(x = 20.7, y = 2.5, xend = 26.5, yend = 2.5)) + annotate(geom="text", size = 5,x=23.5, y= 2.7, label="Extracellular regulators") + scale_y_continuous(limits=c(-3,7.5),breaks=c(-2.5,0,2.5,5.0,7.5))



ggsave("figure_6.pdf",units="in",width= 9,height=7,useDingbats=FALSE)

#figure 7


	#see file transposon-aligning_genes.sh
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

#this is figure 7
ggplot(new_trans_dat_df, aes(x = stage, y = transcript_count)) + geom_sina(aes(colour=transposon),size=0.25,scale="width", alpha=0.25) + stat_summary(aes(group=transposon),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#C75DAB","#009B9E"),name="Transposon-\naligning?") + theme_cowplot() + ylab("Transformed transcript count") + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) + xlab("Stage") + scale_y_continuous(limits=c(-3,20),breaks=c(0,5,10,15,20)) 


ggsave("figure_7.pdf",units="in",width= 9,height=7,useDingbats=FALSE)



#supplemental figure 1, pca


	#see deseq2.R for how this is made
pca_dat <- read.table("all_genes_pca.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)

pca_dat$Stage <- factor(pca_dat$Stage, levels=c("L3", "L4", "adult"))
levels(pca_dat$Stage)[levels(pca_dat$Stage)=="adult"] <- "Adult"

levels(pca_dat$Species)[levels(pca_dat$Species)=="elegans"] <- "C. elegans"
levels(pca_dat$Species)[levels(pca_dat$Species)=="inopinata"] <- "C. inopinata"

ggplot(pca_dat, aes(x = PC1, y = PC2)) + geom_point(aes(colour=Stage,shape=Species),size=2.5) + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12,face="italic"),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("PC1 (49%)") +ylab("PC2 (19%)") + scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73")) 


ggsave("supplemental_figure_1.pdf",units="in",width= 7,height=5.5,useDingbats=FALSE)




#supplemental figure 2


dat <- read.table("rlog_transformed_transcript_counts.tsv", sep="\t", header=TRUE)


head(dat,3)

ino_L3 <- dat[,c(1,2:6)]

ino_L4 <- dat[,c(1,7:10)]

ino_adult <- dat[,c(1,11:15)]

elg_L3 <- dat[,c(1,16:19)]

elg_L4 <- dat[,c(1,20:24)]

elg_adult <- dat[,c(1,25:29)]

ino_L3$ino.L3 <- rowMeans(ino_L3[,2:6], na.rm = TRUE)

ino_L4$ino.L4 <- rowMeans(ino_L4[,2:5], na.rm = TRUE)


ino_adult$ino.adult <- rowMeans(ino_adult[,2:6], na.rm = TRUE)


elg_L3$ele.L3 <- rowMeans(elg_L3[,2:5], na.rm = TRUE)


elg_L4$ele.L4 <- rowMeans(elg_L4[,2:6], na.rm = TRUE)


elg_adult$ele.adult <- rowMeans(elg_adult[,2:6], na.rm = TRUE)


mean_df <- data.frame(gene_id = ino_L3$gen_id,ino_L3 =ino_L3$ino.L3, ino_L4= ino_L4$ino.L4, ino_adult =ino_adult$ino.adult, ele_L3 =elg_L3$ele.L3, ele_L4 = elg_L4$ele.L4, ele_adult = elg_adult$ele.adult)


my_x_title <- expression(paste("Mean ",italic("C. elegans"), " transcript count"))
my_y_title <- expression(paste("Mean ",italic("C. inopinata"), " transcript count "))



a <- ggplot(mean_df, aes(x = ele_L3, y = ino_L3)) + geom_point(size=0.25,alpha=0.25, colour="#E69F00") + theme_cowplot() + theme(axis.title.x=element_blank()) + ylab(my_y_title) + ggtitle("L3")


b <- ggplot(mean_df, aes(x = ele_L4, y = ino_L4)) + geom_point(size=0.25,alpha=0.25, colour= "#56B4E9") + theme_cowplot() + theme(axis.title.y=element_blank())+ xlab(my_x_title) + ggtitle("L4")

c <- ggplot(mean_df, aes(x = ele_adult, y = ino_adult)) + geom_point(size=0.25,alpha=0.25, colour= "#009E73") + theme_cowplot() + theme(axis.title.y=element_blank(),axis.title.x=element_blank()) + ggtitle("Adult")

a+b+c

ggsave("supplemental_figure_2.pdf",units="in",width= 9,height=7,useDingbats=FALSE)

ggsave("supplemental_figure_2.png",units="in",width= 9,height=7,bg="white")
