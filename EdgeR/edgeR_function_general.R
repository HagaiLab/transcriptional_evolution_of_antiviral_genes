# load libararies
library(limma)
library(edgeR)
library(zoo)
library(lmtest)
#library(fdrtools)
library(reshape2)
library(ggplot2)
library(ggsignif)
# Sep 2020

edgeR_func <- function(species_to_analyze, study, unstim_name, stim_name, primate_full_name, rodent_full_name){
  #reading the count matrix! 
  # species_to_analyze_count_data = 'mus_and_homo' - in the new ones- 'homo_vs_mus'
  #species_to_analyze_metadata = 'homo_mus'
  #species_to_analyze_results_name = 'homo_vs_mus'
  
  
  countsMatrix_path = paste("/Users/TzachiHNB1/Documents/docs_from_lilac/tzachi_hagai/project_Sept5_2020_updated/primates_redents_analysis/count_data_for_two_species/",species_to_analyze,'_countData_',study,'.txt',sep='')
  countsMatrix<-read.table(countsMatrix_path,
                         header=T, sep = "\t",row.names=1)
  
  #reading the infotable
  
  info_table_path = paste("/Users/TzachiHNB1/Documents/docs_from_lilac/tzachi_hagai/project_Sept5_2020_updated/primates_redents_analysis/metaData_for_treatments/small_metaData_",study,"_",species_to_analyze,".txt",sep='')
  info_table<-read.table(info_table_path,
                         header = T, sep='\t' ,row.names = 1)
  
  #treatment - making mini-tables with only samples with the "treatment"
  info_table_stim = info_table[info_table$treatment==stim_name,]
  info_table_unstim = info_table[info_table$treatment==unstim_name,]
  
  info_table_primate = info_table[info_table$scientific_name==primate_full_name,]
  info_table_rodent = info_table[info_table$scientific_name==rodent_full_name,]
  
  
  countsMatrix = countsMatrix[,rownames(info_table)]
  
  # round values as salmon gives expected reads rather than actual read count
  countsMatrix<-round(countsMatrix) 
  
  #prints the head of the dataframe
  head(countsMatrix)
  
  # REMOVE LOWLY EXPRESSED GENES 
  countsMatrix_clean = countsMatrix[rowSums(countsMatrix>0)>3,] # remove genes not expressed in at least 3 conditions
  
  
  #mini countMatrix clean, for each treat 
  countsMatrix_stim_clean = countsMatrix_clean[,rownames(info_table_stim)]
  countsMatrix_unstim_clean = countsMatrix_clean[,rownames(info_table_unstim)]
  
  
  # convert to data frame
  
  #general matrix
  countsMatrix_clean = as.data.frame(countsMatrix_clean)
  
  
  #treatment - general
  countsMatrix_stim_clean = as.data.frame(countsMatrix_stim_clean)
  countsMatrix_unstim_clean = as.data.frame(countsMatrix_unstim_clean)
  
  
  # (C) create a DGE object
  #general
  edgeR_DE_obj <- DGEList(counts=countsMatrix_clean, group=info_table$scientific_name)
  
  #treatment- general
  edgeR_DE_obj_stim <- DGEList(counts=countsMatrix_stim_clean, group=info_table_stim$scientific_name)
  edgeR_DE_obj_unstim <- DGEList(counts=countsMatrix_unstim_clean, group=info_table_unstim$scientific_name)
  
  
  
  
  
  # (E)
  # Apply TMM normalization to account for the composition biases
  #general
  edgeR_DE_obj <- calcNormFactors(edgeR_DE_obj, method="TMM") # calc norm factors to adjust raw library counts 
  
  #treatment- general
  edgeR_DE_obj_stim <- calcNormFactors(edgeR_DE_obj_stim, method="TMM") # calc norm factors to adjust raw library counts 
  edgeR_DE_obj_unstim <- calcNormFactors(edgeR_DE_obj_unstim, method="TMM") # calc norm factors to adjust raw library counts 
  
  
  
  # sanity check of TMM norm - do most genes fall around zero?
  par(mar=rep(2,4))
  plotMD(cpm(edgeR_DE_obj, log=TRUE), column=1)
  abline(h=0, col="red", lty=2, lwd=2)
  
  
  
  #edgeR_DE_obj <- estimateDisp(edgeR_DE_obj)
  
  
  edgeR_DE_obj_stim = estimateDisp(edgeR_DE_obj_stim)
  edgeR_DE_obj_unstim <- estimateDisp(edgeR_DE_obj_unstim)
  
  print('levels:')
  print(levels(edgeR_DE_obj_stim$samples$group))
  print(levels(edgeR_DE_obj_unstim$samples$group))
  
  
  et_stim <- exactTest(edgeR_DE_obj_stim, pair=c(primate_full_name,rodent_full_name))
  et_unstim <- exactTest(edgeR_DE_obj_unstim, pair=c(primate_full_name,rodent_full_name))
  
  print('top tags')
  print(topTags(et_stim))
  print(topTags(et_unstim))
  
  #add FDR (q-value) - correction of pval according to benjamine hochberg correction
  #5hours
  et_stim$table$q_val = p.adjust(et_stim$table[,3], "BH")
  et_unstim$table$q_val = p.adjust(et_unstim$table[,3], "BH")
  
  
  # order table of exact test by qval
  #
  et_stim$table <- et_stim$table[order(et_stim$table$q_val),]
  et_unstim$table <- et_unstim$table[order(et_unstim$table$q_val),]
  
  
  # write output
  
  
  f_stim= paste("/Users/TzachiHNB1/Documents/docs_from_lilac/tzachi_hagai/project_Sept5_2020_updated/primates_redents_analysis/edgeR/results_",study,'/',study,'_',stim_name,'_',species_to_analyze,'.txt',sep='')
  f_unstim= paste("/Users/TzachiHNB1/Documents/docs_from_lilac/tzachi_hagai/project_Sept5_2020_updated/primates_redents_analysis/edgeR/results_",study,'/',study,'_',unstim_name,'_',species_to_analyze,'.txt',sep = '')
  
  print(f_stim)
  write.table(et_stim$table, file=f_stim)
  write.table(et_unstim$table, file=f_unstim)
  
  
}
  
#edgeR_func(species_to_analyze = 'homo_vs_mus', study = 'Hagai',unstim_name = 'LF4', stim_name = 'PIC4', primate_full_name = 'Homo sapiens', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'homo_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Homo sapiens', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'pongo_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Pongo abelii', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'pan paniscus_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Pan paniscus', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'gorilla_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Gorilla gorilla', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'pan troglodytes_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Pan troglodytes', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'saimiri_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Saimiri boliviensis boliviensis', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'macaca mulatta_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Macaca mulatta', rodent_full_name = 'Mus musculus')

edgeR_func(species_to_analyze = 'macaca nemestrina_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Macaca nemestrina', rodent_full_name = 'Mus musculus')

# not done!!!
edgeR_func(species_to_analyze = 'papio anubis_vs_mus', study = 'Ploss',unstim_name = 'mock', stim_name = 'treated', primate_full_name = 'Papio anubis', rodent_full_name = 'Mus musculus')

