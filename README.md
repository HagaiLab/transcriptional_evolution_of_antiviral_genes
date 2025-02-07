Transcriptional evolution of antiviral genes between primates and rodents
============================================================================================
This directory contains data files and scripts related to our work: "Evolutionary divergence of induced versus constitutive antiviral gene expression between primates and rodents".
## Folders' details:
1. Input_files folder contains commonly used input files for most of the analyses in the folders bellow, such as human and mouse orthologous genes, lists of the 4 gene group as defined in the manuscript, etc.
   
2. 4_gene_group_expression_patterns folder contains a script for generating plots representing the four groups of human-mouse transcriptionally divergent genes and their patterns of expression. Input files for the script are in input_files folder.
   
3. EdgeR folder contains input files and scripts for executing RNAseq-based DE genes analysis between species (mouse and human) and/or between conditions (unstimulated and stimulated). Further explanations and instructions are found the README file in the folder itself.

4. ChIP_seq folder contains raw ChIP-seq data and downstream analysis for both human and mouse. Further explanations and instructions are found in the README file in the folder itself.

5. Evolutionary_analysis folder contains scripts for evolutionary characteristics analyses of the 4 gene groups: gene age, average dN/dS ,etc. Input files for the script are in input_files folder.

6. TPM_heatmap folder contains Salmon-computed TPM raw data and a script for plotting a heatmap of gene expression patterns for the 4 gene groups.

7. Primates_DE_data contains all DE analysis results of both cross-species datasets used in this manuscript.
