This folder contains raw ChIP-seq data from both mouse and human and an scripts to reproduce the analysis done in the manuscript.

chip_seq_raw folder contains raw and merged NarrowPeak files for human and mouse, as well as bedtools commands for merging the results.

analysis folder contains scripts for downstram analysis of the ChIP-seq raw results as done in the manuscript. In the analysis we also use bedtools toolset. 

Order of ChIP-seq analysis: 
1. Generate .bed files for promoter regions (creating_GTF_and_BED_files.ipynb)
2. If necessary, merge the peaks of chipseq results (bedtools_command_for_merging_narrowPeak_files.txt)
3. Intersect promoter .bed files with merged peaks files (run_bedtools_intersect_script.py)
4. Summarize the results into an excel file -(summaring_chip_seq_results-for_all_1_1_and_DE.ipynb)
5. Use the summarized results as input for chip_seq_analysis.ipynb and chip_seq_results_figures.ipynb 
