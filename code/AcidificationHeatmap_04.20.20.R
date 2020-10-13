
library(tidyverse)
library(viridis)
library(pheatmap)
library(superheat)
library(readxl)
library(colorspace)

homology_table <- read_delim('data/human_homology/rattlesnake_human.BestBlast_results.wSymbol.txt',delim = ' ',col_names = F) %>% 
  mutate(crovirID = str_replace_all(X1,'rattlesnake[|]','')) %>% 
  select(crovirID,humanID=2,symbol=3)



norm.counts <- read_csv('data/norm_counts/cvv_VenomRNAseq_NormCounts_03.19.19.csv') %>% 
  mutate(crovirID = str_split_fixed(X1,'_',2)[,2]) %>% 
  mutate(crovirID = str_replace_all(crovirID,'transcript','protein')) %>% 
  left_join(homology_table) %>% 
  mutate(mergedID = paste(X1,symbol,sep = ':')) %>% 
  select(1,20,21,22,23,16,17,18,19,2,3,4,5,12,13,14,15,6,7,8,9,10,11) %>% 
  unique()


# candidates <- read_xlsx('data/acidification/GeneList_ProtonTransmemTransportGOterm.xlsx',col_names = F)
candidates <- read_xlsx('data/acidification/GeneList_phReductionGOterm.xlsx',col_names = F)


cand.normCounts <- norm.counts %>%
  filter(symbol %in% candidates$...1) %>% 
  arrange(symbol)


cand.normCounts.heat <- as.data.frame(cand.normCounts[,-c(1,2,3,4,5)])
row.names(cand.normCounts.heat) <- cand.normCounts$mergedID

pheatmap(log10(cand.normCounts.heat+1),
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         scale = 'row',
         border_color = 'NA',
         color=magma(50),
         treeheight_row = 5,
         gaps_col = c(4,12),cutree_rows = 4,
)


unext_vs_nonVen <- read_xlsx('data/pairwise_results/significant/phys_only/excel/Unext_vs_NonVenom_Significant.xlsx') %>% 
  mutate(Human = str_split_fixed(Human,'[.]',2)[,1]) %>% 
  left_join(homology_table,by=c('Human'='humanID')) %>% 
  filter(symbol %in% candidates$...1) %>% 
  filter(log2FoldChange > 0) %>% 
  mutate(mergedID = paste(X1,symbol,sep = ':')) %>% 
  unique()

unext_vs_1DPE <- read_xlsx('data/pairwise_results/significant/phys_only/excel/Unext_vs_1DPE_Significant.xlsx') %>% 
  mutate(Human = str_split_fixed(Human,'[.]',2)[,1]) %>% 
  left_join(homology_table,by=c('Human'='humanID')) %>% 
  filter(symbol %in% candidates$...1) %>% 
  filter(log2FoldChange > 0) %>% 
  mutate(mergedID = paste(X1,symbol,sep = ':')) %>% 
  unique()


cand.normCounts.annot <- as.data.frame(rownames(cand.normCounts.heat))
cand.normCounts.annot$Unext_vs_1DPE <- ifelse(cand.normCounts.annot[,1] %in% unext_vs_1DPE$mergedID,1,0)
cand.normCounts.annot$Unext_vs_NonVen <- ifelse(cand.normCounts.annot[,1] %in% unext_vs_nonVen$mergedID,1,0)
row.names(cand.normCounts.annot) <- cand.normCounts.annot[,1]
cand.normCounts.annot <- cand.normCounts.annot[,-1]

pheatmap(log10(cand.normCounts.heat+1),
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         scale = 'row',
         border_color = 'NA',
         color=magma(50),
         treeheight_row = 0,
         gaps_col = c(4,12),
         # cutree_rows = 2,
         annotation_row = cand.normCounts.annot
)


cand.normCounts.ATP <- cand.normCounts %>% 
  filter(substr(symbol,1,1)=='A')

cand.normCounts.other <- cand.normCounts %>% 
  filter(substr(symbol,1,1)!='A')

cand.normCounts.ATP.heat <- as.data.frame(cand.normCounts.ATP[,-c(1,2,3,4,5)])
row.names(cand.normCounts.ATP.heat) <- cand.normCounts.ATP$mergedID

cand.normCounts.other.heat <- as.data.frame(cand.normCounts.other[,-c(1,2,3,4,5)])
row.names(cand.normCounts.other.heat) <- cand.normCounts.other$mergedID


pheatmap(log10(cand.normCounts.ATP.heat+10),
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         scale = 'row',
         border_color = 'NA',
         color=magma(50),
         treeheight_row = 0,
         gaps_col = c(4,12),
         # cutree_rows = 2,
         annotation_row = cand.normCounts.annot
)


pheatmap(log10(cand.normCounts.other.heat+10),
         cluster_cols = F,
         cluster_rows = T,
         cellwidth = 10,
         cellheight = 8,
         fontsize_row = 6,
         scale = 'row',
         border_color = 'NA',
         color=magma(50),
         treeheight_row = 0,
         gaps_col = c(4,12),
         # cutree_rows = 2,
         annotation_row = cand.normCounts.annot
)



