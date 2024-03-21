#Analyzing species-specific amino acid replacements among Caenorhabditis species.

#For Woodruff et al. "Widespread changes in gene expression accompany body size evolution in nematodes"

#set the working directory
setwd("/Users/gavin/genome/transcriptome_revisions_12-2023/SSAR/species-specific_aa_mac_laptop")

#get the species-specific amino acid replacement (SSAAR) data in there, see species-specific_amino_acid_replacements_notes.sh for how this was generated (https://github.com/gcwoodruff/inopinata_developmental_transcriptomics_2023/blob/main/G3_revisions_1/species-specific_amino_acid_replacements/species-specific_amino_acid_replacements_notes.sh)
aadat <- read.csv("number_of_species-specific_aa_changes_for_single_copy_orthologs_b.csv",header=TRUE,stringsAsFactors=TRUE)

keydat <- read.csv("eleg_gene_id_prot_length_fa_file_id_repaired.csv",header=TRUE,stringsAsFactors=TRUE)

#merge SSAAR counts with the key
alldat <- merge(aadat,keydat)

#get the fraction of the protein defined by SSAAR for each species
alldat$brenneri_fra <- alldat$brenneri/alldat$prot_length
alldat$briggsae_fra <- alldat$briggsae/alldat$prot_length
alldat$doughertyi_fra <- alldat$doughertyi/alldat$prot_length
alldat$elegans_fra <- alldat$elegans/alldat$prot_length
alldat$inopinata_fra <- alldat$inopinata/alldat$prot_length
alldat$latens_fra <- alldat$latens/alldat$prot_length
alldat$nigoni_fra <- alldat$nigoni/alldat$prot_length
alldat$remanei_fra <- alldat$remanei/alldat$prot_length
alldat$sinica_fra <- alldat$sinica/alldat$prot_length
alldat$sp26_fra <- alldat$sp26/alldat$prot_length
alldat$sp40_fra <- alldat$sp40/alldat$prot_length
alldat$tropicalis_fra <- alldat$tropicalis/alldat$prot_length
alldat$wallacei_fra <- alldat$wallacei/alldat$prot_length
alldat$kamaaina_fra <- alldat$kamaaina/alldat$prot_length

#get only the numeric data for getting means across rows
just_spsp_aa <- alldat[2:15]

#get the mean number of SSAR for each protein
alldat$mean_spsp_aa <- rowMeans(just_spsp_aa)

#get the SSAR/SSAR_mean ratio for each protein across species
alldat$brenneri_spsp_over_mean <- alldat$brenneri/alldat$mean_spsp_aa
alldat$briggsae_spsp_over_mean <- alldat$briggsae/alldat$mean_spsp_aa
alldat$doughertyi_spsp_over_mean <- alldat$doughertyi/alldat$mean_spsp_aa
alldat$elegans_spsp_over_mean <- alldat$elegans/alldat$mean_spsp_aa
alldat$inopinata_spsp_over_mean <- alldat$inopinata/alldat$mean_spsp_aa
alldat$latens_spsp_over_mean <- alldat$latens/alldat$mean_spsp_aa
alldat$nigoni_spsp_over_mean <- alldat$nigoni/alldat$mean_spsp_aa
alldat$remanei_spsp_over_mean <- alldat$remanei/alldat$mean_spsp_aa
alldat$sinica_spsp_over_mean <- alldat$sinica/alldat$mean_spsp_aa
alldat$sp26_spsp_over_mean <- alldat$sp26/alldat$mean_spsp_aa
alldat$sp40_spsp_over_mean <- alldat$sp40/alldat$mean_spsp_aa
alldat$tropicalis_spsp_over_mean <- alldat$tropicalis/alldat$mean_spsp_aa
alldat$wallacei_spsp_over_mean <- alldat$wallacei/alldat$mean_spsp_aa
alldat$kamaaina_spsp_over_mean <- alldat$kamaaina/alldat$mean_spsp_aa

#how many TOTAL SSAR for each species?

colSums(aadat[-1])

#  brenneri   briggsae doughertyi    elegans  inopinata     latens     nigoni
#      8648       2316       6345      12221      24691       1804       1720
#   remanei     sinica       sp26       sp40 tropicalis   wallacei   kamaaina
#      1972       5440       3082       4571       7095       4126      21148

#summary stats of total number of SSAR

summary(colSums(aadat[-1]))

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1720    2508    5006    7513    8260   24691


#get proteins 20 aa or longer-- I'm interested in the fraction of the protein dominated by SSAR (a kind of normalized SSAR, if you will). So, let's get a respectable per-protein sample size. Arbitrarily choosing 20 aa here.

twendat <- alldat[alldat$prot_length > 19,]

#get summary stats of protein lengths now

summary(twendat$prot_length)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   22.0   166.0   276.0   343.6   436.0  2853.0

#how many proteins now?

nrow(twendat)
#[1] 2767
	#this is the number reported in the results

#new data frame with just fraction of protein defined by SSAR

fracaa <- data.frame(gene=twendat$Public.Name,brenneri = twendat$brenneri_fra,briggsae = twendat$briggsae_fra,doughertyi = twendat$doughertyi_fra,elegans = twendat$elegans_fra,inopinata = twendat$inopinata_fra,latens = twendat$latens_fra,nigoni = twendat$nigoni_fra,remanei = twendat$remanei_fra,sinica = twendat$sinica_fra,sp26 = twendat$sp26_fra,sp40 = twendat$sp40_fra,tropicalis = twendat$tropicalis_fra,wallacei = twendat$wallacei_fra,kamaaina = twendat$kamaaina_fra)

#load the libraries for making the figure
library(reshape2)
library(ggplot2)
library(cowplot)
library(lemon)
library(ggforce)
#melt the data
fracaa_melt <- melt(fracaa,id.vars="gene")

#get the factor levels right -- want it to be in a phylogenetic order
fracaa_melt$variable <- factor(fracaa_melt$variable, levels=c("sp40","sp26","sinica","nigoni","briggsae","remanei","latens","wallacei","tropicalis","doughertyi","brenneri","elegans","inopinata","kamaaina"))

#these species have names now
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="sp40"] <- "tribulationis"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="sp26"] <- "zanzibari"

#add "C. "


levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="tribulationis"] <- "C. tribulationis"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="zanzibari"] <- "C. zanzibari"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="sinica"] <- "C. sinica"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="nigoni"] <- "C. nigoni"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="briggsae"] <- "C. briggsae"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="remanei"] <- "C. remanei"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="latens"] <- "C. latens"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="brenneri"] <- "C. brenneri"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="wallacei"] <- "C. wallacei"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="tropicalis"] <- "C. tropicalis"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="doughertyi"] <- "C. doughertyi"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="elegans"] <- "C. elegans"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="inopinata"] <- "C. inopinata"
levels(fracaa_melt$variable)[levels(fracaa_melt$variable)=="kamaaina"] <- "C. kamaaina"

ggplot(fracaa_melt, aes(x = variable, y = value)) + geom_sina(scale="width",size=0.25,alpha=0.5) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + xlab("Species") + ylab("Fraction of protein with species-specific replacements") + theme_cowplot() + theme(axis.text.y = element_text(face = "italic"))  + coord_flip()

	#this is supplemental figure 8. the phylogeny was added later, see below.

#what are the average SSAAR fractions?

ssaameans <- aggregate(value ~ variable, FUN=mean, data=fracaa_melt)
ssaameans
#           variable       value
#1  C. tribulationis 0.004930711
#2      C. zanzibari 0.003353032
#3         C. sinica 0.006373671
#4         C. nigoni 0.002054957
#5       C. briggsae 0.002539494
#6        C. remanei 0.002268172
#7         C. latens 0.002754567
#8       C. wallacei 0.004557249
#9     C. tropicalis 0.008453768
#10    C. doughertyi 0.006803027
#11      C. brenneri 0.009166554
#12       C. elegans 0.012619901
#13     C. inopinata 0.027058260
#14      C. kamaaina 0.023155561

#SSAAR fraction summary stats
summary(fracaa_melt$value)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000000 0.000000 0.002899 0.008292 0.009950 0.680000

#SSAAR fraction SD
sd(fracaa_melt$value)
#[1] 0.01649807

#Okay, are there really exceptional inopinata proteins, with a super high SSAAR fraction?
#get three standard deviations of mean
mean(fracaa_melt$value)+(sd(fracaa_melt$value)*3)
#0.05778627

#let's mark those inopinata proteins with a SSAAR fraction >3 SD of mean
twendat$high_inop_spsp <- ifelse(twendat$inopinata_fra > 0.05778627,TRUE,FALSE)

#get those genes and remove NA
inop_high_fra <- na.omit(twendat[twendat$high_inop_spsp == TRUE,])
nrow(inop_high_fra)
#[1] 233

#what fraction of the total protein set is this?
nrow(inop_high_fra)/nrow(twendat)
#[1] 0.08420672

#write the genes to a file
#order by inopinata SSAAR fraction
inop_high_fra_order <- inop_high_fra[order(inop_high_fra$inopinata_fra, decreasing = TRUE), ]

#this is just those 233 genes
write(as.character(inop_high_fra_order$Public.Name), "inopinata_genes_high_species-specific_aa_replacements.txt")

#let's get those genes plus some columns of interest
write.table(inop_high_fra_order[c(20,5,6,17,35,24,25,39,40)],"inopinata_genes_high_species-specific_aa_replacements_I.tsv",quote=FALSE,sep='\t')
	#this is supplemental table 21

#write all of the other results

write.table(twendat,"results_R.tsv",row.names=FALSE,quote=FALSE,sep='\t')
	#not in the supplement, but on github



#a phylogeny for terminal branch lengths and for the supplemental figure

#load some libraries
library(ggtree)
library(ape)

#this is the bayesian tree from Stevens et al. 2019 ; retrieved from (https://zenodo.org/record/1402254)
tree <- read.tree("PhyloBayes_species_tree.nwk")

#get only the tips we need
tree_1 <- keep.tip(tree,c("CKAMA","CELEG","CSP34","CWALL","CTROP","CDOUG","CBREN","CREMA","CLATE","CSINI","CSP40","CSP26","CBRIG","CNIGO"))

#root the tree

tree_2 <- root(tree_1, outgroup= "CKAMA")

#new tip labels

tree_2$tip.label <- c("C. kamaaina","C. inopinata","C. elegans","C. remanei","C. latens", "C. tribulationis", "C. zanzibari", "C. sinica","C. nigoni","C. briggsae","C. wallacei","C. tropicalis","C. doughertyi","C. brenneri")


#rotate the tree such that it lines up with our figure
tree_3 <- ggtree(tree_2,branch.length=0.75) + geom_tiplab() + geom_treescale() + xlim(0,1)  
tree_4 <- ggtree::rotate(tree_3,15)
tree_5 <- ggtree::rotate(tree_4, 25)
tree_6 <- ggtree::rotate(tree_5, 20)
tree_7 <- ggtree::rotate(tree_6, 24)
tree_8 <- ggtree::rotate(tree_7, 23)
tree_8
	#this is the tree that was used for the supplemental figure.
	#It was saved to a pdf and then it was merged with the sina plot in illustrator for supplemental figure 8.


#get the terminal branch lengths


n<-length(tree_2$tip.label)
ee<-setNames(tree_2$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree_2$edge[,2])],tree_2$tip.label)
bldf <- data.frame(species=names(ee),branch_lengths=ee)

#merge terminal branch lengths with average fraction of SSAAR
ssaameans
names(ssaameans)[names(ssaameans) == "variable"] <- "species"
names(ssaameans)[names(ssaameans) == "value"] <- "ssaar_mean"
ssablrdf <- merge(ssaameans,bldf)

#make a scatterplot
ggplot(ssablrdf, aes(x = branch_lengths, y = ssaar_mean)) + geom_point() + geom_smooth(method="lm",se=FALSE,size=0.5,linetype="dotted") + xlab("Terminal species branch length") + ylab("Mean fraction of species-specific\namino acid replacements") + theme_cowplot() + scale_y_continuous(limits=c(0,0.03),breaks=c(0,0.0075,0.015,0.0225,0.03))

	#this is supplemental figure 9

#get linear model stats

summary(lm(ssaar_mean~branch_lengths,data=ssablrdf))

#Call:
#lm(formula = ssaar_mean ~ branch_lengths, data = ssablrdf)
#
#Residuals:
       #Min         1Q     Median         3Q        Max
#-0.0027671 -0.0012299  0.0006827  0.0011375  0.0034829
#
#Coefficients:
                 #Estimate Std. Error t value Pr(>|t|)
#(Intercept)    -0.0002076  0.0007856  -0.264    0.796
#branch_lengths  0.1664868  0.0116837  14.249 6.98e-09 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.001913 on 12 degrees of freedom
#Multiple R-squared:  0.9442,	Adjusted R-squared:  0.9395
#F-statistic:   203 on 1 and 12 DF,  p-value: 6.98e-09

#and that's it for SSAAR!!

