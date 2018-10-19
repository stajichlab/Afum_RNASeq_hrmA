# I think that we can slim down the heat map by using
# only the AF293 and hrmA REV strains. 
# I think it would be best to have the replicates collapsed.
# Also the “bad” replicate for R-EV in hypoxia should be dropped. 
# I would also like the normoxia for each together and the hypoxia for 
# each together. 

# As far as what stats to base the heat map off – I think that it should
# be those that are significantly changed (top 100). Previously the heat map 
# said 100 based on p-value – but can we make this so that the changes are 
# log2 > 2  or log2 < -2. Of all the heat maps I like this one a best (below). 
# Can we make the colors more contrasting as well? Lots of requests, sorry!

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
library(GenomicFeatures)
library(genefilter)

if( file.exists("AfumAf293.FungiDB35.db")){
	txdb = loadDb("AfumAf293.FungiDB35.db")
} else {
txdb <- makeTxDbFromGFF("genome/FungiDB-35_AfumigatusAf293.gff",
                        dataSource="FungiDB-35",
                        organism="Aspergillus fumigatus")
saveDb(txdb,"AfumAf293.FungiDB35.db")
}
ebg <- exonsBy(txdb, by="gene")

gene_lengths = sum(width(reduce(ebg)))

my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

countdata <- read.table("results/gsnap_subread/hrmA.Af293.gsnap_reads.tab", header=TRUE, row.names=1)
colnames(countdata) <- gsub("\\.gsnap_Afum_Af293\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("aln\\.", "", colnames(countdata))

countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)
head(countdata)

# this subset are only Af293 and hrmA_REV and
# removed rep2 of hrmA_REV_Hypoxia
samples <- read.csv("samples_plotsubset1.csv",header=TRUE)
exprnames <- do.call(paste,c(samples[c("Strain","Condition","Replicate")],sep="_"))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)

all(exprnames %in% colnames(countdata))
# reorder the columns
countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

# DEseq2 analyses
geno = factor( samples$Strain)
geno = relevel(geno,"AF293")
rep = factor( samples$Replicate)
rep = relevel(rep, "1")

treatment = factor (samples$Condition)
treatment = relevel(treatment,"Hypoxia")

sampleTable <- data.frame(genotype = geno,
                          replicate = rep,
                          condition = treatment
                          )
rownames(sampleTable) = exprnames

dds.oxy <- DESeqDataSetFromMatrix(countData = countdata,
                              colData   = sampleTable, 
                              design    = ~condition)
nrow(dds.oxy)
dds.oxy <- dds.oxy[ rowSums(counts(dds.oxy)) > 1, ]
nrow(dds.oxy)

mcols(dds.oxy)$basepairs = gene_lengths[rownames(mcols(dds.oxy, use.names=TRUE))]

# USE FPM OR FPKM?
# remove genes where the FPKM is too low - I assume we then don't care for 


dds.oxy <- estimateSizeFactors(dds.oxy)
dds.oxy <- estimateDispersions(dds.oxy)

vsd.oxy <- vst(dds.oxy, blind=FALSE)
rld.oxy <- rlog(dds.oxy, blind=FALSE)
head(assay(vsd.oxy), 3)
head(assay(rld.oxy), 3)

# test
dds.oxy <- DESeq(dds.oxy)
res.oxy <- results(dds.oxy)

FPKM <- fpkm(dds.oxy)

min_profile <- data.frame(apply(FPKM, 1, min) )
rownames(min_profile) <- rownames(FPKM)
colnames(min_profile) <- c("Value")
head(min_profile)
filter_genes <- subset(min_profile,min_profile$Value > 5)
filter_res.oxy <- subset(res.oxy,rownames(res.oxy) %in% rownames(filter_genes))
nrow(filter_res.oxy)

select <- order(filter_res.oxy$padj,decreasing=TRUE)[1:50]

df2 <- as.data.frame(colData(dds.oxy)[,c("condition","genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Condition","Genotype")
pdf("plots/heatmaps_showing_replicates.pdf")
pheatmap(assay(vsd.oxy)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Oxy VSD Top Expression")

pheatmap(assay(rld.oxy)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Oxy RLD Top Expression")

# end test

res.oxy <- results(dds.oxy, contrast=c("condition","Normoxia","Hypoxia"))

resLFC.oxy <- lfcShrink(dds.oxy, coef=2, type="apeglm")
summary(resLFC.oxy)

# Get diff expressed
resLFC.oxy <- subset(resLFC.oxy,rownames(resLFC.oxy) %in% rownames(filter_genes))

resSig <- subset(resLFC.oxy, resLFC.oxy$padj < 0.01 & 
                   abs(resLFC.oxy$log2FoldChange) > 2)
resSig <- resSig[order(resSig$padj,decreasing=FALSE),]
write.csv(resSig,"reports/subset1_condition_Normoxia_vs_Hypoxia.csv")
write.csv(fpm(dds.oxy),"reports/subset1_FPM.csv")
write.csv(fpkm(dds.oxy),"reports/subset1_FPKM.csv")

vsd.oxy.reorder <- assay(vsd.oxy)[,c("AF293_Normoxia.r1","AF293_Normoxia.r2","AF293_Normoxia.r3",
                                     "hrmA_REV_Normoxia.r1","hrmA_REV_Normoxia.r2","hrmA_REV_Normoxia.r3",
                                     "AF293_Hypoxia.r1","AF293_Hypoxia.r2","AF293_Hypoxia.r3",
                                     "hrmA_REV_Hypoxia.r1","hrmA_REV_Hypoxia.r3")]

rld.oxy.reorder <- assay(rld.oxy)[,c("AF293_Normoxia.r1","AF293_Normoxia.r2","AF293_Normoxia.r3",
                                     "hrmA_REV_Normoxia.r1","hrmA_REV_Normoxia.r2","hrmA_REV_Normoxia.r3",
                                     "AF293_Hypoxia.r1","AF293_Hypoxia.r2","AF293_Hypoxia.r3",
                                     "hrmA_REV_Hypoxia.r1","hrmA_REV_Hypoxia.r3")]

pheatmap(rld.oxy.reorder[rownames(resSig),], cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD p<0.01 and 2-fold Exp Fold change Oxy")

pheatmap(vsd.oxy.reorder[rownames(resSig),], cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD p<0.01 and 2-fold Exp Fold change Oxy")

# == now test for difference when using genotype
dds.geno <- DESeqDataSetFromMatrix(countData = countdata,
                                   colData   = sampleTable, 
                                   design    = ~genotype)

nrow(dds.geno)
dds.geno <- dds.geno[ rowSums(counts(dds.geno)) > 1, ]
nrow(dds.geno)

mcols(dds.geno)$basepairs = gene_lengths[rownames(mcols(dds.geno, use.names=TRUE))]

dds.geno <- estimateSizeFactors(dds.geno)
dds.geno <- estimateDispersions(dds.geno)

vsd.geno <- vst(dds.geno, blind=FALSE)
rld.geno <- rlog(dds.geno, blind=FALSE)
head(assay(vsd.geno), 3)
head(assay(rld.geno), 3)

# test
dds.geno <- DESeq(dds.geno)
res.geno <- results(dds.geno)

filter_res.geno <- subset(res.geno,rownames(res.geno) %in% rownames(filter_genes))

select <- order(filter_res.geno$padj,
                decreasing=TRUE)

df2 <- as.data.frame(colData(dds.geno)[,c("condition","genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Condition","Genotype")
pheatmap(assay(vsd.geno)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Geno VSD Top Expression Diff")
pheatmap(assay(rld.geno)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Geno RLD Top Expression Diff")

# end test

res.geno <- results(dds.geno, contrast=c("genotype","AF293","hrmA_REV"))

resLFC.geno <- lfcShrink(dds.geno, coef=2, type="apeglm")
summary(resLFC.geno)

# Get diff expressed
resLFC.geno <- subset(resLFC.geno,rownames(resLFC.geno) %in% rownames(filter_genes))

resSig <- subset(resLFC.geno, resLFC.geno$padj < 0.05 & 
                   abs(resLFC.geno$log2FoldChange) > 2)
resSig <- resSig[order(resSig$padj,decreasing=FALSE),]
write.csv(resSig,"reports/subset1_genotype_AF293_vs_hrmA_REV.csv")
# not needed- these are same no matter the experimental design
# write.csv(fpm(dds.geno),"reports/subset1_FPM_Genotype.csv")
# write.csv(fpkm(dds.geno),"reports/subset1_FPKM_Genotype.csv")

vsd.geno.reorder <- assay(vsd.geno)[,c("AF293_Normoxia.r1","AF293_Normoxia.r2","AF293_Normoxia.r3",
                                     "hrmA_REV_Normoxia.r1","hrmA_REV_Normoxia.r2","hrmA_REV_Normoxia.r3",
                                     "AF293_Hypoxia.r1","AF293_Hypoxia.r2","AF293_Hypoxia.r3",
                                     "hrmA_REV_Hypoxia.r1","hrmA_REV_Hypoxia.r3")]

rld.geno.reorder <- assay(rld.geno)[,c("AF293_Normoxia.r1","AF293_Normoxia.r2","AF293_Normoxia.r3",
                                     "hrmA_REV_Normoxia.r1","hrmA_REV_Normoxia.r2","hrmA_REV_Normoxia.r3",
                                     "AF293_Hypoxia.r1","AF293_Hypoxia.r2","AF293_Hypoxia.r3",
                                     "hrmA_REV_Hypoxia.r1","hrmA_REV_Hypoxia.r3")]

pheatmap(rld.geno.reorder[rownames(resSig),], cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD p<0.01 and 2-fold Exp Fold change Geno")

pheatmap(vsd.geno.reorder[rownames(resSig),], cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD p<0.01 and 2-fold Exp Fold change Geno")



# == now generate collapsed
dds.oxy$id <- factor(paste0(dds.oxy$genotype,'.',dds.oxy$condition),
                      levels = c("AF293.Normoxia",
                                 "hrmA_REV.Normoxia",
                                 "AF293.Hypoxia",
                                 "hrmA_REV.Hypoxia"))


dds.oxy.Coll <- collapseReplicates(dds.oxy,dds.oxy$id)
colData(dds.oxy.Coll)


dds.oxy.Coll <- DESeq(dds.oxy.Coll)
res.oxy.Coll <- results(dds.oxy.Coll)
rld.oxy.Coll <- rlog(dds.oxy.Coll, blind=FALSE)

mcols(dds.oxy.Coll)$basepairs = gene_lengths[rownames(mcols(dds.oxy.Coll, use.names=TRUE))]

resSig.oxy.Coll <- subset(resSig.oxy.Coll,rownames(resSig.oxy.Coll) %in% rownames(filter_genes))
resSig.oxy.Coll <- subset(res.oxy.Coll,res.oxy$padj < 0.01)

summary(resSig.oxy.Coll)

resSig.oxy.Coll <- resSig.oxy.Coll[order(resSig.oxy.Coll$pvalue,decreasing=FALSE),]

write.csv(resSig.oxy.Coll,"reports/subset1_Collapsed_Normox_vs_Hypoxia.csv")
write.csv(fpm(dds.oxy.Coll),"reports/subset1_Collapsed_FPM.csv")
write.csv(fpkm(dds.oxy.Coll),"reports/subset1_Collapsed_FPKM.csv")

resBest.oxy.Coll <- subset(resSig.oxy.Coll,abs(resSig.oxy.Coll$log2FoldChange) > 2)

mat.oxy.Coll <- rld.oxy.Coll[ rownames(resBest.oxy.Coll), ]
colconditions.oxy.Coll = as.data.frame(colData(dds.oxy.Coll)[,c("condition","genotype")])

#mat.Coll.reorder <- assay(mat.Coll)[,c("AF293.Normoxia","hrmA_REV.Normoxia",
#                                     "AF293.Hypoxia","hrmA_REV.Hypoxia")]

pdf("plots/heatmaps_showing_replicates.pdf",width=10)
pheatmap(assay(mat.oxy.Coll),method="complete",
         main = "Collapsed Reps - Oxygen model, p-value < 0.01 and log_fold_change > 2", 
         show_colnames=TRUE, show_rownames = FALSE,
         annotation_legend = TRUE, legend=TRUE, cluster_rows=TRUE, 
         cluster_cols = FALSE, cexRow=0.3,
         annotation_col=colconditions.oxy.Coll)

resBest.oxy.Coll4 <- subset(resSig.oxy.Coll,abs(resSig.oxy.Coll$log2FoldChange) > 4)

mat.oxy.Coll4 <- rld.oxy.Coll[ rownames(resBest.oxy.Coll4), ]
colconditions.oxy.Coll = as.data.frame(colData(dds.oxy.Coll)[,c("condition","genotype")])

pheatmap(assay(mat.oxy.Coll4),method="complete",
         main = "Collapsed Reps - Oxygen model, p-value < 0.01 and log_fold_change > 4", 
         show_colnames=TRUE, show_rownames = FALSE,
         annotation_legend = TRUE, legend=TRUE, cluster_rows=TRUE, 
         cluster_cols = FALSE, cexRow=0.3,
         annotation_col=colconditions.oxy.Coll)


resBest.oxy.Coll4 <- subset(resSig.oxy.Coll,abs(resSig.oxy.Coll$log2FoldChange) > 4)
pheatmap(assay(mat.oxy.Coll4),method="complete",
         main = "Collapsed Reps - Oxygen model, p-value < 0.01 and log_fold_change > 4 top 100", 
         show_colnames=TRUE, show_rownames = TRUE,
         annotation_legend = TRUE, legend=TRUE, cluster_rows=TRUE, 
         cluster_cols = FALSE, cexRow=0.3,
         annotation_col=colconditions.oxy.Coll)
# == dds geno
dds.geno$id <- factor(paste0(dds.geno$genotype,'.',dds.geno$condition),
                     levels = c("AF293.Normoxia",
                                "hrmA_REV.Normoxia",
                                "AF293.Hypoxia",
                                "hrmA_REV.Hypoxia"))


dds.geno.Coll <- collapseReplicates(dds.geno,dds.geno$id)
colData(dds.geno.Coll)

dds.geno.Coll <- DESeq(dds.geno.Coll)
res.geno.Coll <- results(dds.geno.Coll)
rld.geno.Coll <- rlog(dds.geno.Coll, blind=FALSE)

mcols(dds.geno.Coll)$basepairs = gene_lengths[rownames(mcols(dds.geno.Coll, use.names=TRUE))]

resSig.geno.Coll <- subset(res.geno.Coll,res.geno$padj < 0.01)
resSig.geno.Coll <- resSig.geno.Coll[order(resSig.geno.Coll$padj,decreasing=FALSE),]

write.csv(resSig.geno.Coll,"reports/subset1_Collapsed_AF293_vs_hrmA_REV.csv")

resBest.geno.Coll <- subset(resSig.geno.Coll,abs(resSig.geno.Coll$log2FoldChange) > 2)

mat.geno.Coll <- rld.geno.Coll[ rownames(resBest.geno.Coll), ]
colconditions.geno.Coll = as.data.frame(colData(dds.geno.Coll)[,c("condition","genotype")])

#mat.Coll.reorder <- assay(mat.Coll)[,c("AF293.Normoxia","hrmA_REV.Normoxia",
#                                     "AF293.Hypoxia","hrmA_REV.Hypoxia")]
pheatmap(assay(mat.geno.Coll),method="complete",
         main = "Collapsed Reps - Geno model, p-value < 0.01 and log_fold_change > 2", 
         show_colnames=TRUE, show_rownames = TRUE,
         annotation_legend = TRUE, legend=TRUE, cluster_rows=TRUE, 
         cluster_cols = FALSE, cexRow=0.3,
         annotation_col=colconditions.geno.Coll)


pdf("plots/PCA_2types_expression.pdf")
pcaData <- plotPCA(vsd.oxy, intgroup=c("genotype","condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# FIX ME HERE
ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=treatment,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=replicate,label=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()

ggplot(pcaData, aes(PC1, PC2, color=genotype,shape=replicate,label=treatment)) +
  geom_point(size=3) + geom_text(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()
dev.off()
