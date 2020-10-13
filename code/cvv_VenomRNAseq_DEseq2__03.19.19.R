
library(DESeq2)
library(pheatmap)
library(viridis)
library(IHW)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)


tx2gene <- read.csv('data/misc/tx2gene_v3.csv',header = F,stringsAsFactors = F)

rawCounts <- read.csv('data/raw_counts/cvv_venomRNAseq_rawCounts_03.19.19.csv',stringsAsFactors = F,skip=1,header=T)
rawCounts <- rawCounts[,c(-2,-3,-4,-5)]

rawCounts.simple <- rawCounts[which(!is.na(rawCounts$Crovir_Transcript_ID)),]

rawCounts.simple <- rawCounts.simple[,c(-1,-2)]

rawCounts.simple <- merge(rawCounts.simple,tx2gene,by=1,all.x=T)

row.names(rawCounts.simple) <- rawCounts.simple$V2
rawCounts.simple <- rawCounts.simple[,c(-1,-20)]

#rawCounts.simple <- rawCounts.simple[which(rowSums(rawCounts.simple) > 1),]
rawCounts.simple <- rawCounts.simple[rowSums( rawCounts.simple != 0 ) > 8,]


colData <- (DataFrame(condition=group <- factor(c('1DPE','1DPE','1DPE','1DPE','NonVen','NonVen','NonVen','NonVen','NonVen','NonVen','3DPE','3DPE','3DPE','3DPE','Unext','Unext','Unext','Unext'))))


dds <- DESeqDataSetFromMatrix(rawCounts.simple,colData,formula(~condition))

dds <- DESeq(dds)

normcounts <- counts(dds,normalized=TRUE)
#write.csv(as.data.frame(normcounts),file='./cvv_VenomRNAseq_NormCounts_03.19.19.csv',row.names = T)
maxCooks <- as.data.frame(mcols(dds)$maxCooks)
allcooks <- assays(dds)[["cooks"]]

hist(maxCooks[,1],breaks = 10000,xlim=c(0,6))

###
### Variance stabilized PCA
###

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=3) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_minimal()


###
### Unextracted vs. 1DPE
###

deResUnextVs1DPE <- as.data.frame(results(dds, cooksCutoff = 10, contrast=c('condition','1DPE','Unext')))

ihwRes <- ihw(pvalue ~ baseMean,  data = deResUnextVs1DPE, alpha = 0.05)

rejections(ihwRes)
deRes <- deResUnextVs1DPE
deRes <- na.omit(deRes)

# plot(ihwRes)
# 
# gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
#   geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
# gg
# gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)
# 
# ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
# 
# deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 10)
# 
# ggplot(deRes, aes(x=pvalue)) +
#   geom_histogram(binwidth = 0.025, boundary = 0) +
#   facet_wrap( ~ baseMeanGroup, nrow = 2)
# 
# ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")
# 
# rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
#                  covariate_type="base mean"),
#       data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
#                  covariate_type="log2 fc")) %>%
#   ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
#   facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))


deResUnextVs1DPE$IHW_pvalue <- ihwRes@df$adj_pvalue

deResUnextVs1DPE <- deResUnextVs1DPE[order(deResUnextVs1DPE$IHW_pvalue),]

# write.csv(as.data.frame(deResUnextVs1DPE),file='./pairwise_results/cvv_Unextracted.vs.1DPE_VG_PairwiseResult_03.28.19.csv',row.names = T)


###
### 1DPE vs. 3DPE
###

deRes1DPEvs3DPE <- as.data.frame(results(dds, cooksCutoff = 10, contrast=c('condition','3DPE','1DPE')))

ihwRes <- ihw(pvalue ~ baseMean,  data = deRes1DPEvs3DPE, alpha = 0.05)

rejections(ihwRes)
deRes <- deRes1DPEvs3DPE
deRes <- na.omit(deRes)


deRes1DPEvs3DPE$IHW_pvalue <- ihwRes@df$adj_pvalue

deRes1DPEvs3DPE <- deRes1DPEvs3DPE[order(deRes1DPEvs3DPE$IHW_pvalue),]

# write.csv(as.data.frame(deRes1DPEvs3DPE),file='./pairwise_results/cvv_1DPE.vs.3DPE_VG_PairwiseResult_03.28.19.csv',row.names = T)



###
### Unext vs. NonVenom
###

deResUnextVsNonVen <- as.data.frame(results(dds, cooksCutoff = 10, contrast=c('condition','Unext','NonVen')))

ihwRes <- ihw(pvalue ~ baseMean,  data = deResUnextVsNonVen, alpha = 0.05)

rejections(ihwRes)
deRes <- deResUnextVsNonVen
deRes <- na.omit(deRes)

deResUnextVsNonVen$IHW_pvalue <- ihwRes@df$adj_pvalue

deResUnextVsNonVen <- deResUnextVsNonVen[order(deResUnextVsNonVen$IHW_pvalue),]

# write.csv(as.data.frame(deResUnextVsNonVen),file='./pairwise_results/cvv_Unext.vs.NonVenom_VG_PairwiseResult_03.29.19.csv',row.names = T)


        