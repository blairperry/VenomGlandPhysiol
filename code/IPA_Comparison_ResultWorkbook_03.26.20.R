
library(tidyverse)
library(colorspace)
library(readxl)
library(janitor)




# Canonical Pathways ------------------------------------------------------

# cpa.z <- read_xls('data/ipa_results/Comparison/AllSig_CPAcomparison_zscore.xls',skip = 1, na = 'N/A') %>% 
#   select(pathway=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) %>% 
#   pivot_longer(-1, names_to = 'analysis', values_to = 'zscore')
# 
# 
# cpa.p <- read_xls('data/ipa_results/Comparison/AllSig_CPAcomparison_pval.xls',skip = 1, na = 'N/A') %>% 
#   select(pathway=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) %>% 
#   pivot_longer(-1, names_to = 'analysis', values_to = 'pval')
# 
# cpa.full <- cpa.z %>% left_join(cpa.p)
# 
# 
# cpa.full %>% 
#   filter(pval > 1.3) %>% 
#   group_by(pathway) %>% 
#   mutate(totalZ = sum(zscore)) %>% 
#   filter(min(zscore) > 0 & totalZ >= 1) %>% 
#   ggplot(aes(x=reorder(pathway,totalZ),y=zscore,fill=analysis)) +
#   geom_bar(stat='identity',show.legend = T,color='white') +
#   scale_fill_manual(values = c('#8F8F8F','seagreen')) +
#   # facet_wrap(~analysis) +
#   ylab('Activation Z-Score') +
#   xlab('') +
#   coord_flip() +
#   # theme_bw() +
#   theme_classic() +
#   theme(panel.grid = element_blank())
# 

# v2

cpa.z2 <- read_xls('data/ipa_results/Comparison/AllSig_CPAcomparison_zscore.xls',skip = 1, na = 'N/A') %>% 
  select(pathway=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) 


cpa.p2 <- read_xls('data/ipa_results/Comparison/AllSig_CPAcomparison_pval.xls',skip = 1, na = 'N/A') %>% 
  select(pathway=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) 

cpa.full2.paths <- cpa.z2 %>% 
  left_join(cpa.p2, by = c("pathway"), suffix = c('.z','.p')) %>% 
  filter(Unext_vs_NonVen.z >= 1 & Unext_vs_NonVen.p > 1.3) %>% 
  select(1)

cpa.z2 <- cpa.z2 %>% 
  pivot_longer(-1, names_to = 'analysis', values_to = 'zscore')
cpa.p2 <- cpa.p2 %>% 
  pivot_longer(-1, names_to = 'analysis', values_to = 'pval')

cpa.full2 <- cpa.full2.paths %>% 
  left_join(cpa.z2) %>% 
  left_join(cpa.p2) %>% 
  filter(pval > 1.3)

cpa.full2 %>% 
  group_by(pathway) %>% 
  mutate(totalZ = sum(zscore)) %>% 
  # filter(min(zscore) > 0 & totalZ >= 1) %>% 
  ggplot(aes(x=reorder(pathway,totalZ),y=zscore,fill=analysis)) +
  geom_bar(stat='identity',show.legend = T,color='white') +
  scale_fill_manual(values = c('#8F8F8F','seagreen')) +
  # facet_wrap(~analysis) +
  ylab('Activation Z-Score') +
  xlab('') +
  coord_flip() +
  # theme_bw() +
  theme_classic() +
  theme(panel.grid = element_blank())

# URMs ------------------------------------------------------

# urm.z <- read_xls('data/ipa_results/Comparison/AllSig_URMcomparison_zscore.xls',skip = 1, na = 'N/A') %>% 
#   select(urm=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) %>% 
#   pivot_longer(-1, names_to = 'analysis', values_to = 'zscore')
# 
# 
# urm.p <- read_xls('data/ipa_results/Comparison/AllSig_URMcomparison_pval.xls',skip = 1, na = 'N/A') %>% 
#   select(urm=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) %>% 
#   pivot_longer(-1, names_to = 'analysis', values_to = 'pval')
# 
# urm.full <- urm.z %>% left_join(urm.p)
# 
# urm.full %>% 
#   filter(pval > 1.3) %>% 
#   group_by(urm) %>% 
#   mutate(totalZ = sum(zscore)) %>% 
#   filter(min(zscore) > 0 & totalZ >= 3) %>%
#   ggplot(aes(x=reorder(urm,totalZ),y=zscore,fill=analysis)) +
#   geom_bar(stat='identity',show.legend = T,color='white') +
#   scale_fill_manual(values = c('salmon1','skyblue1')) +
#   # facet_wrap(~analysis) +
#   ylab('Activation Z-Score') +
#   xlab('') +
#   coord_flip() +
#   # theme_bw() +
#   theme_classic() +
#   theme(panel.grid = element_blank())

## v2

urm.z2 <- read_xls('data/ipa_results/Comparison/AllSig_URMcomparison_zscore.xls',skip = 1, na = 'N/A') %>% 
  select(urm=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) 


urm.p2 <- read_xls('data/ipa_results/Comparison/AllSig_URMcomparison_pval.xls',skip = 1, na = 'N/A') %>% 
  select(urm=1,Unext_vs_NonVen=2,Unext_vs_1DPE=3) 

urm.full2.paths <- urm.z2 %>% 
  left_join(urm.p2, by = c("urm"), suffix = c('.z','.p')) %>% 
  filter(Unext_vs_NonVen.p > 1.3) %>% 
  # top_n(25,wt = Unext_vs_NonVen.z) %>%
  filter(Unext_vs_NonVen.z >= 2) %>% 
  select(1)

urm.z2 <- urm.z2 %>% 
  pivot_longer(-1, names_to = 'analysis', values_to = 'zscore')
urm.p2 <- urm.p2 %>% 
  pivot_longer(-1, names_to = 'analysis', values_to = 'pval')

urm.full2 <- urm.full2.paths %>% 
  left_join(urm.z2) %>% 
  left_join(urm.p2) %>% 
  filter(pval > 1.3)

urm.full2 %>% 
  mutate(urm = str_split_fixed(urm,'[ (]',2)[,1]) %>% 
  group_by(urm) %>% 
  mutate(totalZ = sum(zscore)) %>% 
  # filter(min(zscore) > 0 & totalZ >= 1) %>% 
  ggplot(aes(x=reorder(urm,totalZ),y=zscore,fill=analysis)) +
  geom_bar(stat='identity',show.legend = T,color='white') +
  scale_fill_manual(values = c('#8F8F8F','seagreen')) +
  # facet_wrap(~analysis) +
  ylab('Activation Z-Score') +
  xlab('') +
  coord_flip() +
  # theme_bw() +
  theme_classic() +
  theme(panel.grid = element_blank())

