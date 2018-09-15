library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(fdrtool)
library(geneplotter)
library(EDASeq)

my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

countdata <- read.table("reports/Hypoxia.Af293.gsnap_reads.nostrand.tab", header=TRUE, row.names=1)
colnames(countdata) <- gsub("\\.gsnap_Afum_Af293\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("aln\\.", "", colnames(countdata))
countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
head(countdata)

samples <- read.csv("samples.csv",header=TRUE)
exprnames <- do.call(paste,c(samples[c("Strain","Condition","Replicate")],sep="_"))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)

# check that experimental columns match in order
all(exprnames %in% colnames(countdata))
all(exprnames == colnames(countdata))

# reorder the columns anyways... in case data change along the way
countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

# DEseq2 analyses
geno = factor( samples$Strain)
rep = factor( samples$Replicate)
treatment = factor (samples$Condition)

sampleTable <- data.frame(genotype = geno,
                          replicate = rep,
                          condition = treatment)
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData   = sampleTable, 
                              design    = genotype ~ condition )

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/Hypoxia_RNASeq_gsnap.pdf")

plotDispEsts(dds)

multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

MA.idx = t(combn(1:4, 2))
for( i in  seq_along( MA.idx[,1])){ 
  MDPlot(counts(dds, normalized = T), 
         c(MA.idx[i,1],MA.idx[i,2]), 
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
}

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]


#datCollapsed <- collapseReplicates(dds, groupby=dds$genotype,run=dds$replicate,renameCols=TRUE)

df2 <- as.data.frame(colData(dds)[,c("genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Genotype")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")

topVar <- head(order(rowVars(assay(vsd)),
          decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
        fontsize_row = 7,fontsize_col = 7,
        cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")

pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Top Expression")

topVar <- head(order(rowVars(assay(rld)),
    decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Most different")



sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=treatment)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() + theme_bw()


norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

topVarGenes <- order(-rowVars(log.norm.counts)[0:1000])
mat<-log.norm.counts[topVarGenes,]
mat<-mat -rowMeans(mat)

pheatmap(mat,method="complete",main = "TopVar normalized", show_rownames = F,
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.55 )

topVarGenes1 <- order(-rowVars(assay(rld)))[0:1000]
mat1 <- assay(rld)[ topVarGenes1, ]
mat1<- mat1 - rowMeans(mat1)

pheatmap(mat1, method="complete",
         main = "Unsupervised 1000 genes ",
         show_rownames = F,annotation_legend = FALSE, 
         legend=T, cluster_cols=TRUE)
           
dds <- DESeq(dds)
res <- results(dds)
res    
res <- results(dds, contrast=c("condition","Normoxia","Hypoxia"))
res

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_Normoxia_vs_Hypoxia", type="apeglm")
resLFC
summary(resLFC)
summary(res)
res05 <- results(dds,alpha=0.05)
summary(res05)
sum(res$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))


