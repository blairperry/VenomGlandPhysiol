
library(tidyverse)

##
## Gene Sets to Generate for GO Analysis
##

unext.vs.nonven <- read_csv('data/pairwise_results/significant/phys_only/Unext_vs_NonVenom_Significant.csv')

unext.vs.1DPE <- read_csv('data/pairwise_results/significant/phys_only/Unext_vs_1DPE_Significant.csv')

o1DPE.vs.3DPE <- read_csv('data/pairwise_results/significant/phys_only/1DPE_vs_3DPE_Significant.csv')



### Background (all expressed genes that went into DEseq2)

tx2gene <- read.csv('data/misc/tx2gene_v4.csv',header = F,stringsAsFactors = F)
rawCounts <- read.table('data/misc/cvv_venomRNAseq_rawCounts_updatedGFF_11.19.19.txt',stringsAsFactors = F,skip=1,header=T)
rawCounts <- rawCounts[,c(-2,-3,-4,-5)]

rawCounts.simple <- rawCounts[which(!is.na(rawCounts$Crovir_Transcript_ID)),]
rawCounts.simple <- rawCounts.simple[,c(-1,-2)]
row.names(rawCounts.simple) <- rawCounts.simple$Crovir_Transcript_ID
rawCounts.simple <- rawCounts.simple[,-1]
rawCounts.simple <- rawCounts.simple[rowSums( rawCounts.simple != 0 ) > 8,]

allGeneBackground <- rawCounts.simple %>% 
  rownames_to_column(var = 'Crovir_ID') %>% 
  select(Crovir_ID) %>% 
  mutate(Crovir_ProtID = str_replace_all(Crovir_ID,'transcript','protein'))

homology_table <- read_delim('data/human_homology/rattlesnake_human.BestBlast_results.txt',delim='\t',col_names = F) %>% 
  mutate(Crovir_ProtID = str_split_fixed(X1,'[|]',2)[,2]) %>% 
  mutate(X2 = str_split_fixed(X2,'[|]',2)[,2]) %>% 
  select(Crovir_ProtID,HumanID = X2)


allGeneBackground.humID <- allGeneBackground %>% 
  left_join(homology_table)

#write.csv(allGeneBackground.humID,'data/pairwise_results/significant/phys_only/_subsets/Background_forGOanalysis.csv',row.names = F)


### Unext vs. Non-Venom (Unique)

postExt <- union(unext.vs.1DPE$X1,o1DPE.vs.3DPE$X1)

unext.vs.nonven_unique <- unext.vs.nonven %>% 
  filter(!(X1 %in% postExt))

#write.csv(unext.vs.nonven_unique,'data/pairwise_results/significant/phys_only/_subsets/Unext_vs_NonVenom_UNIQUE_Sig.csv',row.names = F)


### Unext vs. Non-Venom (Unique Upreg in VG)

unext.vs.nonven_unique_up <- unext.vs.nonven_unique %>% 
  filter(log2FoldChange > 0)
  
#write.csv(unext.vs.nonven_unique_up,'data/pairwise_results/significant/phys_only/_subsets/Unext_vs_NonVenom_UP_UNIQUE_Sig.csv',row.names = F)


### Unext vs. Non-Venom (Unique Downreg in VG)

unext.vs.nonven_unique_dwn <- unext.vs.nonven_unique %>% 
  filter(log2FoldChange < 0)

#write.csv(unext.vs.nonven_unique_dwn,'data/pairwise_results/significant/phys_only/_subsets/Unext_vs_NonVenom_DWN_UNIQUE_Sig.csv',row.names = F)


### Overlap: Unext vs. Non-Venom Upreg and  Unext vs 1DPE Upreg 

overlap <- unext.vs.nonven %>% 
  filter(X1 %in% unext.vs.1DPE$X1) %>% 
  filter(!(X1 %in% o1DPE.vs.3DPE$X1)) %>% 
  left_join(unext.vs.1DPE,suffix = c('_Unext.vs.NonVen','_Unext.vs.1DPE'),by='X1')

#write.csv(overlap,'data/pairwise_results/significant/phys_only/_subsets/Unext_vs_1DPE_vs_NonVen_OVERLAP.csv',row.names = F)


### Unext vs. 1DPE (Unique)

filterout <- union(unext.vs.nonven$X1,o1DPE.vs.3DPE$X1)

unext.vs.1DPE_unique <- unext.vs.1DPE %>% 
  filter(!(X1 %in% filterout))

#write.csv(unext.vs.1DPE_unique,'data/pairwise_results/significant/phys_only/_subsets/Unext_vs_1DPE_UNIQUE.csv',row.names = F)



### Unext vs. 1DPE (Unique Upreg)

unext.vs.1DPE_unique_up <- unext.vs.1DPE_unique %>% 
  filter(log2FoldChange > 0)

#write.csv(unext.vs.1DPE_unique_up,'data/pairwise_results/significant/phys_only/_subsets/Unext_vs_1DPE_UP_UNIQUE.csv',row.names = F)


### Unext vs. 1DPE (Unique Downreg)

unext.vs.1DPE_unique_dwn <- unext.vs.1DPE_unique %>% 
  filter(log2FoldChange < 0)


#write.csv(unext.vs.1DPE_unique_dwn,'data/pairwise_results/significant/phys_only/_subsets/Unext_vs_1DPE_DWN_UNIQUE.csv',row.names = F)
