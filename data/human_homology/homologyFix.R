
library(tidyverse)

homology_table <- read_delim('data/human_homology/rattlesnake_human.BestBlast_results.txt',delim='\t',col_names = F) %>% 
  mutate(X2 = str_replace_all(X2,'human[|]','')) %>% 
  mutate(X2 = str_split_fixed(X2,'[.]',2)[,1])

human_annot <- read_csv('data/human_homology/proteins_51_582967.csv') %>% 
  select('Accession'=9,'symbol'=7) %>% 
  mutate(Accession = str_split_fixed(Accession,'[.]',2)[,1])


homology_table2 <- homology_table %>% 
  left_join(human_annot,by=c('X2'='Accession')) %>% 
  mutate(symbol = ifelse(is.na(symbol),X2,symbol))


write.table(homology_table2,'data/human_homology/rattlesnake_human.BestBlast_results.wSymbol.txt',col.names = F,row.names = F,quote = F)
