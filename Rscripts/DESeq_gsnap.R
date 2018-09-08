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

countdata <- read.table("Af293_Hypoxia.gsnap_reads.tab", header=TRUE, row.names=1)
colnames(countdata) <- gsub("\\.gsnap_Afum_Af293\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("^aln/", "", colnames(countdata))

samples <- read.csv("samples.csv",header=TRUE)
exprnames <- do.call(paste,c(samples[c("Strain","Condition","Replicate")],sep="."))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)

# DEseq2 analyses
geno = factor( samples$Strain)
rep = factor( samples$Replicate)
treatment = factor (samples$Condition)

sampleTable <- data.frame(genotype = geno,
                          condition = treatment
                          rep = replicate)
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromMatrix(coundata,sampleTable, treatment ~ genotype)

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/Hypoxia_RNASeq_kallisto.pdf")

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


datCollapsed <- collapseReplicates(dds, groupby=dds$genotype,run=dds$replicate,renameCols=TRUE)

df2 <- as.data.frame(colData(dds)[,c("genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Genotype")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD")

topVar <- head(order(rowVars(assay(vsd)),
          decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]
mat  <- mat - rowMeans(mat)

mat2  <- assay(vsd)[ topVar, ]

controlAve <- rowMeans(mat2[ , Genotype == "AF293" ])

pheatmap(mat, show_rownames=TRUE,
        fontsize_row = 7,fontsize_col = 7,
        cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")

pheatmap(datCollapsed, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Normalized to Af293")




pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD")

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

ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
