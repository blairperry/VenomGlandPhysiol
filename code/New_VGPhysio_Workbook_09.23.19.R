
library(tidyverse)
library(cowplot)





# Loading human homologs ------------------------------------------------

homolog_table <- read_delim('data/human_homology/rattlesnake_human.BestBlast_results.txt',delim = '\t',col_names = F) %>% 
  mutate(Rattlesnake = str_split_fixed(X1,'\\|',2)[,2]) %>% 
  mutate(Rattlesnake = str_replace(Rattlesnake,'protein','transcript')) %>% 
  mutate(Human = str_split_fixed(X2,'\\|',2)[,2]) %>% 
  select(Rattlesnake,Human)


# Read in all pairwise result files ---------------------------------------------

read_wName <- function(flnm) {
  read_csv(flnm,na = 'NA') %>% 
    mutate(test = str_split_fixed(flnm,'_',4)[[3]])
}


allData <-
  list.files(path='./data/pairwise_results',pattern='*.csv',full.names = T,recursive = F) %>% 
  map_df(~read_wName(.)) %>% 
  mutate(test = str_remove_all(test,'\\.')) %>% 
  mutate(t_ID = str_split_fixed(X1,'_',2)[,2]) %>% 
  left_join(homolog_table,by=c('t_ID' = 'Rattlesnake'))


signif <- allData %>% 
  filter(IHW_pvalue < 0.05)

unique(signif$test)

sum(is.na(allData$Human)) / length(allData$Human)


Unext_1DPE <- signif %>% 
  filter(test == 'Unextractedvs1DPE') %>% 
  mutate(type = ifelse(str_detect(X1,'Venom'),'Venom','Nonvenom')) %>% 
  select(-9)

Unext_3DPE <- signif %>% 
  filter(test == 'Unextvs3DPE') %>% 
  mutate(type = ifelse(str_detect(X1,'Venom'),'Venom','Nonvenom')) %>% 
  select(-9)

OneDPE_3DPE <-  signif %>% 
  filter(test == '1DPEvs3DPE') %>% 
  mutate(type = ifelse(str_detect(X1,'Venom'),'Venom','Nonvenom')) %>% 
  select(-9)

OneDPE_NonVenom <-  signif %>% 
  filter(test == '1DPEvsNonVenom') %>% 
  mutate(type = ifelse(str_detect(X1,'Venom'),'Venom','Nonvenom')) %>% 
  select(-9)

ThreeDPE_NonVenom <-  signif %>% 
  filter(test == '3DPEvsNonVenom') %>% 
  mutate(type = ifelse(str_detect(X1,'Venom'),'Venom','Nonvenom')) %>% 
  select(-9)

Venom_NonVenom <-  signif %>% 
  filter(test == 'VenomvsNonVenom') %>% 
  mutate(type = ifelse(str_detect(X1,'Venom'),'Venom','Nonvenom')) %>% 
  select(-9)

Unext_NonVenom <-  signif %>% 
  filter(test == 'UnextvsNonVenom') %>% 
  mutate(type = ifelse(str_detect(X1,'Venom'),'Venom','Nonvenom')) %>% 
  select(-9)

#write_csv(Unext_1DPE,'./data/pairwise_results/significant/Unext_vs_1DPE_Significant.csv')
#write_csv(Unext_3DPE,'./data/pairwise_results/significant/Unext_vs_3DPE_Significant.csv')
#write_csv(OneDPE_3DPE,'./data/pairwise_results/significant/1DPE_vs_3DPE_Significant.csv')
#write_csv(OneDPE_NonVenom,'./data/pairwise_results/significant/1DPE_vs_NonVenom_Significant.csv')
#write_csv(ThreeDPE_NonVenom,'./data/pairwise_results/significant/3DPE_vs_NonVenom_Significant.csv')
#write_csv(Venom_NonVenom,'./data/pairwise_results/significant/Venom_vs_NonVenom_Significant.csv')
#write_csv(Unext_NonVenom,'./data/pairwise_results/significant/Unext_vs_NonVenom_Significant.csv')



# Filter out venom genes to make physio gene sets -------------------------

phys_Unext_1DPE <- Unext_1DPE %>% 
  filter(type == 'Nonvenom')

phys_Unext_3DPE <- Unext_3DPE %>% 
  filter(type == 'Nonvenom')

phys_OneDPE_3DPE <- OneDPE_3DPE %>% 
  filter(type == 'Nonvenom')

phys_OneDPE_NonVenom <- OneDPE_NonVenom %>% 
  filter(type == 'Nonvenom')

phys_ThreeDPE_NonVenom <- ThreeDPE_NonVenom %>% 
  filter(type == 'Nonvenom')

phys_Venom_NonVenom <- Venom_NonVenom %>% 
  filter(type == 'Nonvenom')

phys_Unext_NonVenom <- Unext_NonVenom %>% 
  filter(type == 'Nonvenom')



# write_csv(phys_Unext_1DPE,'./data/pairwise_results/significant/phys_only/Unext_vs_1DPE_Significant.csv')
# write_csv(phys_Unext_3DPE,'./data/pairwise_results/significant/phys_only/Unext_vs_3DPE_Significant.csv')
# write_csv(phys_OneDPE_3DPE,'./data/pairwise_results/significant/phys_only/1DPE_vs_3DPE_Significant.csv')
# write_csv(phys_OneDPE_NonVenom,'./data/pairwise_results/significant/phys_only/1DPE_vs_NonVenom_Significant.csv')
# write_csv(phys_ThreeDPE_NonVenom,'./data/pairwise_results/significant/phys_only/3DPE_vs_NonVenom_Significant.csv')
# write_csv(phys_Venom_NonVenom,'./data/pairwise_results/significant/phys_only/Venom_vs_NonVenom_Significant.csv')
#write_csv(phys_Unext_NonVenom,'./data/pairwise_results/significant/phys_only/Unext_vs_NonVenom_Significant.csv')



