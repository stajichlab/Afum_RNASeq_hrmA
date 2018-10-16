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

countdata <- read.table("results/gsnap_subread/hrmA.Af293.gsnap_reads.tab", header=TRUE, row.names=1)
colnames(countdata) <- gsub("\\.gsnap_Afum_Af293\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("aln\\.", "", colnames(countdata))

countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
head(countdata)

samples <- read.csv("samples.csv",header=TRUE)
samples = samples[-c(23),] # remove a bad sample
exprnames <- do.call(paste,c(samples[c("Strain","Condition","Replicate")],sep="_"))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)

#exprnames = exprnames[!exprnames %in% "hrmA_REV_Hypoxia.r2"]

# check that experimental columns match in order
all(exprnames %in% colnames(countdata))
#all(exprnames == colnames(countdata))
# reorder the columns
countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

# DEseq2 analyses
geno = factor( samples$Strain)
rep = factor( samples$Replicate)
treatment = factor (samples$Condition)

sampleTable <- data.frame(genotype = geno,
                          replicate = rep,
                          condition = treatment
)
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData   = sampleTable, 
                                  design    = ~condition)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

pdf("plots/PCA_expression.pdf")
pcaData <- plotPCA(vsd, intgroup=c("genotype","condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# FIX ME HERE
ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=treatment,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw() + ggtitle("VSD PCA plot - no Bad Reps")

ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=replicate,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw() + ggtitle("VSD PCA plot  - no Bad Reps")


ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=replicate,label=treatment)) +
  geom_point(size=3) + geom_text(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw() + ggtitle("VSD PCA plot - no Bad Reps")


