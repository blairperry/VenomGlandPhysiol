
library(tidyverse)
library(readxl)


homology <- read_delim('data/human_homology/rattlesnake_human.BestBlast_results.wSymbol.txt',delim=' ',col_names = F) %>% 
  select('crovir'=1,'Human'=2,'symbol'=3) %>% 
  mutate(crovir = str_split_fixed(crovir,'[|]',2)[,2]) %>% 
  mutate(crovir = str_replace_all(crovir,'protein','transcript'))

unext_vs_nonven <- read_xlsx('data/pairwise_results/significant/phys_only/excel/Unext_vs_NonVenom_Significant.xlsx') %>% 
  filter(log2FoldChange > 0) %>% 
  mutate(Human  = str_split_fixed(Human,'[.]',2)[,1]) %>% 
  left_join(homology) %>% 
  unique()


unext_vs_1dpe <- read_xlsx('data/pairwise_results/significant/phys_only/excel/Unext_vs_1DPE_Significant.xlsx') %>% 
  filter(log2FoldChange > 0) %>% 
  mutate(Human  = str_split_fixed(Human,'[.]',2)[,1]) %>% 
  left_join(homology) %>% 
  unique()

#write.table(unext_vs_nonven,'data/goAnalysis/Unext_vs_NonVen_InputUpregGenes.txt',sep = '\t')
#write.table(unext_vs_1dpe,'data/goAnalysis/Unext_vs_1dpe_InputUpregGenes.txt',sep = '\t')


bgList <- read_csv('data/norm_counts/cvv_VenomRNAseq_NormCounts_03.19.19.csv') %>% 
  mutate(crovir = str_split_fixed(X1,'_',2)[,2]) %>% 
  left_join(homology) %>% 
  unique() %>% 
  select(symbol)

write.table(bgList,'data/goAnalysis/_UPDATED_geneSymbol_bglist.txt',sep = '\t',row.names = F,quote = F)


