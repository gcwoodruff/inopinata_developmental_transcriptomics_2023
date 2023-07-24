

sed '1d' L3_stats_results.tsv > L3_stats_results.tsv
sed '1d' L4_stats_results.tsv > L4_stats_results.tsv
sed '1d' adult_stats_results.tsv > adult_stats_results.tsv

cd all_stages_LRT_results

awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"L3"}' L3_stats_results.tsv > L3_stats_results.tsv.tmp
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"L4"}' L4_stats_results.tsv > L4_stats_results.tsv.tmp
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0,"Adult"}' adult_stats_results.tsv > adult_stats_results.tsv.tmp

cat L3_stats_results.tsv.tmp L4_stats_results.tsv.tmp adult_stats_results.tsv.tmp > all_stages_LRT_results.tsv.tmp


echo -e "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tstage" | cat - all_stages_LRT_results.tsv.tmp > all_stages_LRT_results.tsv

rm L3_stats_results.tsv.tmp
rm L3_stats_results.tsv
rm L4_stats_results.tsv.tmp
rm L4_stats_results.tsv
rm adult_stats_results.tsv.tmp
rm adult_stats_results.tsv
rm all_stages_LRT_results.tsv.tmp

	#gumienny_savage-dunn_tgf-beta_list.txt generated from table 1 of http://wormbook.org/chapters/www_tgfbsignal.2/TGFbetasignal.html
grep -w -f gumienny_savage-dunn_tgf-beta_list.txt all_stages_LRT_results.tsv > all_stages_LRT_results_TGF-beta.tsv.tmp

echo -e "gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tstage" | cat - all_stages_LRT_results_TGF-beta.tsv.tmp > all_stages_LRT_results_TGF-beta.tsv

rm all_stages_LRT_results_TGF-beta.tsv.tmp