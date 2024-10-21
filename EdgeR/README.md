This folder contains the input data and scripts required for reproducing bulk-rna DE analysis performed by EdgeR package.

1. count_matrix contains a count matrix (genes X samples) with TPM counts from the bulk-rna-seq.
2. info_tables contains metadata for each sample (species, unstimulated or stimulated, etc.). Using the info table we can decide what is the base of comparison between the samples.
3. scripts:
     a. edgeR_function_general - an exmaple of running DE analysis between human and mouse within unstimulated or stimulated samples.
     b. egdeR_glm_interactions - an exmaple of running DE analysis using a linear model containing the     
Species:Trearment interaction trem. 
