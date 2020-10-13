
library(janitor)
library(tidyverse)
library(colorspace)
library(readxl)
library(pheatmap)
library(tidytext)


urm_types <- c('transcription regulator','growth factor','kinase','transporter','cytokine','complex','enzyme','G-protein coupled receptor','group','ion channel','ligand-dependent nuclear receptor','peptidase','phosphotase','transmembrane receptor')

# Upstream Regulatory Molecules  ------------------------------------------

urm.unextVs1DPE <- read_excel('data/ipa_results/VenomProduction/new_UnextVs1DPE_AllSig_URMresults_03.25.20.xls', skip = 1, na = 'NA') %>% 
  select(urm=1,expLogRatio=2,type=3,predState=4,zscore=5,pval=6,molecules=7) %>% 
  mutate(zscore = as.numeric(zscore)) %>% 
  filter(pval < 0.05) %>% 
  filter(type %in% urm_types)

top25.urm.unextVs1DPE <- urm.unextVs1DPE %>% 
  # filter(zscore > 0) %>% 
  # arrange(-zscore) %>% 
  # top_n(25,wt=zscore)
  filter(zscore >= 2) %>% 
  arrange(-zscore)


 
ggplot(top25.urm.unextVs1DPE,aes(x=reorder(urm,zscore),y=zscore)) +
  geom_bar(stat='identity',show.legend = F) +
  # scale_fill_continuous_sequential('Mint',name='-Log10(p-value)') +
  # geom_point(aes(x=reorder(urm,zscore),y=-0.25,color=type),size=3.5) +
  # scale_color_discrete_qualitative(palette = 'Set 3',name='URM Type') +
  ylab('Activation Z-Score') +
  ggtitle('Unextracted vs. 1DPE',subtitle = 'Top 25 Upregulated URMs') +
  scale_y_continuous(expand = c(0.01,0)) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        plot.title.position = 'plot'
  )




# Canonical Pathways ------------------------------------

## Unextracted vs. 1DPE

cpa.UnextVs1DPE <- read_excel('data/ipa_results/VenomProduction/new_UnextVs1DPE_AllSig_CPAresults_03.25.20.xls',skip=1,na = 'NaN') %>% 
  mutate(mol_num = str_count(Molecules,',') + 1) %>% 
  mutate(pvalue = 10^(-1 * `-log(p-value)`)) %>% 
  select(Pathway=1,7,zscore=4,3,5,6) %>% 
  filter(pvalue < 0.05) %>% 
  mutate(direction = ifelse(is.na(zscore), 'Undetermined', ifelse(zscore > 0, 'Activated',ifelse(zscore < 0, 'Inhibited','Undetermined'))))

cpa.UnextVs1DPE %>% 
  filter(direction == 'Activated' & zscore >= 1) %>% 
  top_n(20,wt=abs(zscore)) %>% 
  ggplot(aes(x=reorder(Pathway,abs(zscore)),y=abs(zscore))) +
  geom_bar(stat='identity') +
  scale_fill_continuous_sequential('Mint',name='-Log10(p-value)') +
  ggtitle('Unextracted vs. 1DPE Venom Gland') +
  # facet_wrap(~direction,scales = 'free',ncol=1) +
  ylab('Activation Z-Score') +
  scale_y_continuous(expand = c(0.01,0)) +
  coord_flip() +
  theme_classic(base_size = 12) + 
  theme(panel.grid = element_blank(),
        plot.title.position = 'plot',
        axis.title.y = element_blank(),
        strip.background = element_rect(colour = NA),
        legend.position = 'bottom')

