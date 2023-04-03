library("DESeq2")
library(ggplot2)
library(htmltools)
library(dplyr)
library(ggpubr)
library(NMF)

# read the data; counts, descriptor, and annotation
countData <- read.csv('data.csv', header = T, sep =",")
metadata <- read.csv('metadata.csv', header=T, sep=',')
annotation = read.csv('annotation.csv')

# run DESeq2 pipeline
dseq.data = DESeqDataSetFromMatrix(countData=countData, colData=metadata, 
                                   design=~condition, tidy=T )
dseq.data = DESeq(dseq.data)

# join to annotation data, and filter to only significant results
results = as_tibble(results(dseq.data, tidy=T))
results = results %>% arrange(pvalue)
results = results %>% inner_join(annotation, by=c('row'='ensgene'))
sig.results = results %>% filter(padj<0.05)

# plot the counts of the four most significant genes from the control to
# treated group
plot1 = plotCounts(dseq.data, gene=as.character(sig.results[1,1]), intgroup="condition", 
                   returnData=TRUE) %>% 
  ggplot(aes(condition, count)) + 
  geom_boxplot(aes(fill=condition)) + 
  scale_y_log10() + 
  ggtitle(as.character(sig.results[1,1]))
plot2 = plotCounts(dseq.data, gene=as.character(sig.results[2,1]), intgroup="condition", 
                   returnData=TRUE) %>% 
  ggplot(aes(condition, count)) + 
  geom_boxplot(aes(fill=condition)) + 
  scale_y_log10() + 
  ggtitle(as.character(sig.results[2,1]))
plot3 = plotCounts(dseq.data, gene=as.character(sig.results[3,1]), intgroup="condition", 
                   returnData=TRUE) %>% 
  ggplot(aes(condition, count)) + 
  geom_boxplot(aes(fill=condition)) + 
  scale_y_log10() + 
  ggtitle(as.character(sig.results[3,1]))
plot4 = plotCounts(dseq.data, gene=as.character(sig.results[4,1]), intgroup="condition", 
                   returnData=TRUE) %>% 
  ggplot(aes(condition, count)) + 
  geom_boxplot(aes(fill=condition)) + 
  scale_y_log10() + 
  ggtitle(as.character(sig.results[4,1]))

ggarrange(plot1, plot2, plot3, plot4 + rremove("x.text"), labels=c('A', 'B', 'C','D'), nrow=2, ncol=2)

# create volcano plot
results = results %>% mutate(sig=padj<0.05)
ggplot(data=results, aes(x=log2FoldChange, y=-log10(pvalue), col=sig)) +
  geom_point()

# data transform for a heatmap of top 50 DEGs
transformed = vst(dseq.data, blind=F)
aheatmap(assay(transformed)[arrange(results, padj, pvalue)$row[1:50],], labRow=arrange(results, padj, pvalue)$symbol[1:50], scale='row', distfun='pearson', annCol=select(metadata, condition), col=c('green', 'black', 'black', 'red'))

