summary chipseq analysis (ordered): 
1. generate bed files for promoters (creating_promoters_files)
2. merge the peaks of chipseq results or atac (example in commands_in_cluster_fr_merge_narrowPeak)
generated the merged peaks files - human_narrowPeak_files or mouse
3. intersect promoter file with merged_peaks files (example in run_bedtools_intersect_script)
generated the results from cluster
4. summarize the results from cluster - summaring_chip_seq_results-for_all_1_1_and_DE 
this goes to the metadata 
5. do the chip_seq analysis in the tzachi_hagai\thesis\scripts_and_figures\chip_seq and the analysing_summary_chip_seq_only_DE_genes_check script
