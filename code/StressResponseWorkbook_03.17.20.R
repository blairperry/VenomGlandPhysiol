
library(tidyverse)
library(pheatmap)
library(readxl)
library(janitor)
library(colorspace)



# Unfolded Protein Response -----------------------------------------------

## Gene heatmap

upr.fc <- read_xls('data/stress_response/UPR_AllSigComp_foldChange.xls',skip = 1,na = 'N/A') %>%
  clean_names() %>% 
  select(gene=1,Unext_vs_NonVenom=2,Unext_vs_1DPE=3,Unext_vs_3DPE=4) %>% 
  arrange(-Unext_vs_NonVenom)
  
upr.pv <- read_xls('data/stress_response/UPR_AllSigComp_pval.xls',skip = 1) %>% clean_names()         # Currently don't really need since only significant genes went into analysis..

upr.fc.heat <- as.data.frame(upr.fc[,2:4])
row.names(upr.fc.heat) <- upr.fc$gene

paletteLength <- 40
mycols <- diverging_hcl(paletteLength, palette = "Blue-Red 2")

lim <- 5

myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))


pheatmap(upr.fc.heat,
         cluster_cols = F,
         cluster_rows = F,
         treeheight_row = 0,
         cellheight = 12,
         cellwidth = 12,
         border_color = 'white',
         color = mycols,
         breaks = myBreaks,
         main='Genes',
         na_col = 'grey95',
         filename = 'figures/raw/ipa/UPR_GeneHeatmap_03.17.20.pdf')


## URM Heatmap

urm_types <- c('transcription regulator','growth factor','kinase','transporter','cytokine','complex','enzyme','G-protein coupled receptor','group','ion channel','ligand-dependent nuclear receptor','peptidase','phosphotase','transmembrane receptor')

upr.mols <- read_delim('data/stress_response/UPR_allMolecules.txt',delim='\t',col_names = F) %>% filter(str_detect(1,' ',negate=T))

all.urms <- read_xls('data/ipa_results/z_old/urm/allSig/AllSig_Comparison_URM_Table.xls',skip=1) %>% 
# all.urms <- read_xls('data/ipa_results/urm/allSig/AllSig_Comparison_URM_Table.xls',skip=1) %>% 
  clean_names() %>% 
  filter(molecule_type %in% urm_types) %>% 
  select(1,2,8,9)

upr.urms <- all.urms %>% filter(upstream_regulator %in% upr.mols$X1) %>% 
  filter(p_value_of_overlap < 0.05) %>% 
  select(-4) %>% 
  pivot_wider(-1,names_from = analysis,values_from = activation_z_score) %>% 
  select(1,Unext_vs_NonVen = 2, Unext_vs_1DPE = 3, DPE1_vs_3DPE = 4) %>% 
  arrange(-Unext_vs_NonVen)

upr.urms.heat <- as.data.frame(upr.urms[,-1])
row.names(upr.urms.heat) <- upr.urms$upstream_regulator

paletteLength <- 40

mycols <- colorRampPalette(c('#136364','#60C2AF','grey95','#EAA155','#D9442E'))(paletteLength)


lim <- 7

myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

pheatmap(upr.urms.heat,
         cluster_cols = F,
         cluster_rows = F,
         treeheight_row = 10,
         cellheight = 12,
         cellwidth = 12,
         border_color = 'white',
         color = mycols,
         breaks = myBreaks,
         main='URMs',
         na_col = 'grey95',
         filename = 'figures/raw/ipa/UPR_UrmHeatmap_diffColor_03.31.20.pdf'
         )



# Endoplasmic Reticulum Stress Pathway -----------------------------------------------

## Gene heatmap

ersp.fc <- read_xls('data/stress_response/ERSP_AllSigComp_foldchange.xls',skip = 1,na = 'N/A') %>%
  clean_names() %>% 
  select(gene=1,Unext_vs_NonVenom=2,Unext_vs_1DPE=3,Unext_vs_3DPE=4) %>% 
  arrange(-Unext_vs_NonVenom)

ersp.pv <- read_xls('data/stress_response/ERSP_AllSigComp_pval.xls',skip = 1) %>% clean_names()         # Currently don't really need since only significant genes went into analysis..

ersp.fc.heat <- as.data.frame(ersp.fc[,2:4])
row.names(ersp.fc.heat) <- ersp.fc$gene

paletteLength <- 40
mycols <- diverging_hcl(paletteLength, palette = "Blue-Red 2")

lim <- 5

myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))


pheatmap(ersp.fc.heat,
         cluster_cols = F,
         cluster_rows = F,
         cellheight = 12,
         cellwidth = 12,
         border_color = 'white',
         color = mycols,
         breaks = myBreaks,
         main='Genes',
         na_col = 'grey95',
         filename = 'figures/raw/ipa/ERSP_GeneHeatmap_03.17.20.pdf')


## URM Heatmap

ersp.mols <- read_delim('data/stress_response/ERSP_allMolecules',delim='\t',col_names = F) %>% filter(str_detect(1,' ',negate=T))


ersp.urms <- all.urms %>% filter(upstream_regulator %in% ersp.mols$X1) %>% 
  filter(p_value_of_overlap < 0.05) %>% 
  select(-4) %>% 
  pivot_wider(-1,names_from = analysis,values_from = activation_z_score) %>% 
  select(1,Unext_vs_NonVen = 2, Unext_vs_1DPE = 3, DPE1_vs_3DPE = 4) %>% 
  arrange(-Unext_vs_NonVen)

ersp.urms.heat <- as.data.frame(ersp.urms[,-1])
row.names(ersp.urms.heat) <- ersp.urms$upstream_regulator

paletteLength <- 40

mycols <- diverging_hcl(paletteLength, palette = "Blue-Red 2")

lim <- 7

myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

pheatmap(ersp.urms.heat,
         cluster_cols = F,
         cluster_rows = F,
         treeheight_row = 10,
         cellheight = 12,
         cellwidth = 12,
         border_color = 'white',
         color = mycols,
         breaks = myBreaks,
         main='URMs',
         na_col = 'grey95',
         filename = 'figures/raw/ipa/ERSP_UrmHeatmap_03.17.20.pdf')




# Combined pathway heatmaps -----------------------------------------------

## Genes


combGene <- upr.fc %>% 
  bind_rows(ersp.fc) %>% 
  unique() %>% 
  mutate(UPR = ifelse(gene %in% upr.fc$gene,1,0)) %>% 
  mutate(ERSP = ifelse(gene %in% ersp.fc$gene,1,0)) %>% 
  arrange(-Unext_vs_NonVenom)

combGene.heat <- as.data.frame(combGene[,2:4])
row.names(combGene.heat) <- combGene$gene


combGene.annot <- as.data.frame(combGene[,5:6])
row.names(combGene.annot) <- combGene$gene

paletteLength <- 40

mycols <- diverging_hcl(paletteLength, palette = "Blue-Red 2")

lim <- 7

myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

pheatmap(combGene.heat,
         cluster_cols = F,
         cluster_rows = F,
         treeheight_row = 10,
         cellheight = 12,
         cellwidth = 12,
         border_color = 'white',
         color = mycols,
         breaks = myBreaks,
         main='URMs',
         na_col = 'grey95',
         annotation_row = combGene.annot,
         filename = 'figures/raw/ipa/stress_response/combGeneURMHeatmap_03.17.20.pdf')

## Urms

combined <- upr.urms %>% 
  bind_rows(ersp.urms) %>% 
  unique() %>% 
  mutate(UPR = ifelse(upstream_regulator %in% upr.urms$upstream_regulator,1,0)) %>% 
  mutate(ERSP = ifelse(upstream_regulator %in% ersp.urms$upstream_regulator,1,0))

combined.heat <- as.data.frame(combined[,2:4])
row.names(combined.heat) <- combined$upstream_regulator


combined.annot <- as.data.frame(combined[,5:6])
row.names(combined.annot) <- combined$upstream_regulator

paletteLength <- 40

mycols <- diverging_hcl(paletteLength, palette = "Blue-Red 2")

lim <- 7

myBreaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

pheatmap(combined.heat,
         cluster_cols = F,
         cluster_rows = F,
         treeheight_row = 10,
         cellheight = 12,
         cellwidth = 12,
         border_color = 'white',
         color = mycols,
         breaks = myBreaks,
         main='URMs',
         na_col = 'grey95',
         annotation_row = combined.annot,
         filename = 'figures/raw/ipa/stress_response/combinedURMHeatmap_03.17.20.pdf')


