

library(tidyverse)
library(ggrepel)


normCounts <- read_csv('data/norm_counts/cvv_VenomRNAseq_NormCounts_03.19.19.csv') %>% 
  mutate(ID = str_split_fixed(X1,'-',3)[,3])

venList <- read_csv('data/venom_genes/venom_1DPEAvgNormCounts.csv') %>% 
  select(ID = X1)

venCounts <- venList %>% left_join(normCounts) %>% 
  filter(!is.na(.$X1)) %>% 
  filter(ID != 'PLA2_gIIE') %>% 
  pivot_longer(3:20,names_to = 'replicate',values_to = 'normCounts') %>% 
  select(-2) %>% 
  mutate(treatment = ifelse(str_detect(replicate,'ODPE'),'ODPE',
                            ifelse(str_detect(replicate,'TDPE'),'TDPE',
                                   ifelse(str_detect(replicate,'Unext'),'Unext','Body')))) %>% 
  mutate(Family = ifelse(str_detect(ID,'PLA2'),'PLA2',
                            ifelse(str_detect(ID,'SVMP'),'SVMP',
                                   ifelse(str_detect(ID,'SVSP'),'SVSP','Other')))) %>% 
  mutate(normCounts = normCounts+1) %>% 
  mutate(log10normCounts = log10(normCounts)) %>% 
  mutate(high = ifelse(normCounts > 100000,'Yes','No'))

## Defining "Venom relevant venom genes" as those w/ norm counts > 100,000 in at least one time point


venCounts %>% 
  mutate(treatment = factor(treatment, levels=c('Body','Unext','ODPE','TDPE'))) %>% 
  group_by(ID) %>% 
  filter('Yes' %in% high) %>% 
  ungroup() %>% 
  group_by(ID,treatment,Family) %>% 
  summarize(mean = mean(normCounts),stdev = sd(normCounts)) %>% 
  ggplot(aes(x=treatment,y=mean,group=ID,color=Family)) + 
  geom_line() +
  geom_text_repel(aes(label=ID)) +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), width=.05) + 
  geom_point() +
  #scale_y_log10() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())


