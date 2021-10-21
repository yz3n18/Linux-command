
# set working path and creat directory
data.path <- "./"; setwd(data.path) #create dir
res.path    <- file.path(data.path, "Results")
fig.path    <- file.path(data.path, "Figures")

if (!file.exists(data.path)) { dir.create(data.path) }
if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }

library(ggplot2)
library(dplyr)
library(tibble)
library(janitor)

exprSet<-read.table('/Users/yihuawang/Mark_Pax/all.id.txt',sep = '\t',header = T) %>%
  select(1,7:35) %>% # remove W sample
  column_to_rownames('Geneid') 
colnames(exprSet)<-gsub('.bam','',colnames(exprSet))
group_list<-read.xlsx('/Users/yihuawang/Mark_Pax/NovogeneCPETstudyWallis.xlsx')%>%
  slice( -(1:7))%>%
  row_to_names(row_number = 1)
exprSet<-exprSet[,group_list$`Sample Name`]
group_list$Group<-c(rep('IPF_0',11),
                    rep('IPF_8',8),
                    rep('Control',10))
condition_rna<-factor(c(rep('IPF_0',11),
                        rep('IPF_8',8),
                        rep('Control',10)))## make condition
condition_rna


coldata <- data.frame(row.names=colnames(exprSet), condition_rna)
coldata
dds <- DESeqDataSetFromMatrix(countData=exprSet, colData=coldata, design=~condition_rna)### create a matrix for further analysis
nrow(dds)
dds <- DESeq(dds)
res <- results(dds)

png("RAWvsNORM.png")
rld <- rlogTransformation(dds)
exprSet_new=assay(rld)
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition_rna))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(exprSet_new)))


pdf("qc-heatmap-samples.pdf", width =10, height =10)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition_rna], RowSideColors=mycols[condition_rna],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()


exprSet[, c(1:8)] <- sapply(exprSet[, c(1:8)], as.integer) # convert into integer

pca_data <- plotPCA(rld, intgroup=c('condition_rna'), returnData=T, ntop=5000)
?plotPCA
percentVar <- round(100 * attr(pca_data, "percentVar"))
pdf('PCA.pdf',width = 8,height = 8)
ggplot(pca_data,aes(PC1, PC2, color=condition_rna,label=name))+
  geom_point(size=3) +
  #geom_text_repel(data=pca_data,nudge_x = 0.5)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

