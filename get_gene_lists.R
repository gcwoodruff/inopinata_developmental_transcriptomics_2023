
dat <- read.table("time_course_LRT_results_Speciesinopinata.StageL4_by_int_term.tsv", sep="\t", header=TRUE)

dat$int_term_x_neg_one <- dat$log2FoldChange*(-1)

dat$neg.log.p.adj <- -(log10(dat$padj))

dat$int_term_abs <- abs(dat$log2FoldChange)

dat$is.sig <- dat$padj < 0.05

dat.om <- na.omit(dat)


dat.om$is.sig <- as.factor(dat.om$is.sig)
levels(dat.om$is.sig)[levels(dat.om$is.sig)=="FALSE"] <- "No"
levels(dat.om$is.sig)[levels(dat.om$is.sig)=="TRUE"] <- "Yes"

#get lists for enrichment analyses, top 10% of significant genes ordered in various ways



dat$int.direction <- ifelse(dat$int_term_x_neg_one >0, "positive","negative")

dat_sig <- subset(dat, padj < 0.05)

nrow(dat_sig)
#7204

#order by absolute value of interaction term (high to low)

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




#order by p ; same as before



dat_order_int_p <- dat_sig[order(-dat_sig$padj),]

dat_order_int_p_list <- head(dat_order_int_p,720)

write.table(dat_order_int_p_list, file= 'LIST_D_top_ten_perc_p_value_Species--L4-adult_int_term_only_sig.tsv',sep="\t", quote = FALSE)


#top ten abs val

dat_order_abs_val_int_h_l_top_ten <- head(dat_order_abs_val_int_h_l,10)


ggplot(dat_order_abs_val_int_h_l_top_ten, aes(y=reorder(gene, int_term_abs), x=int_term_abs)) + geom_bar(stat='identity',fill="#9ecae1") + theme_cowplot() +xlab("|Interaction Term|") + ylab("Gene")
