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

#txdb <- makeTxDbFromGFF("genome/FungiDB-35_AfumigatusAf293.gff",
#                        dataSource="FungiDB-35",
#                        organism="Aspergillus fumigatus")
#saveDb(txdb,"AfumAf293.FungiDB35.db")
txdb = loadDb("AfumAf293.FungiDB35.db")
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
#all(exprnames == colnames(countdata))
# reorder the columns
countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

# DEseq2 analyses
geno = factor( samples$Strain)
rep = factor( samples$Replicate)
treatment = factor (samples$Condition)

sampleTable <- data.frame(condition = treatment,
                          genotype = geno,
                          replicate = rep
                          )
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData   = sampleTable, 
                              design    = genotype ~ condition )

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

mcols(dds)$basepairs = gene_lengths[rownames(mcols(dds, use.names=TRUE))]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

# test
dds <- DESeq(dds)
res <- results(dds)
select <- order(res$padj,
                decreasing=TRUE)[1:50]
df2 <- as.data.frame(colData(dds)[,c("condition","genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Genotype","Hypoxia")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD Top Expression")
# end test

dds <- DESeq(dds)
res <- results(dds)
res <- results(dds, contrast=c("condition","Normoxia","Hypoxia"))

resLFC <- lfcShrink(dds, coef=2, type="apeglm")
summary(resLFC)

# Get diff expressed
resSig <- subset(resLFC, resLFC$padj < 0.05)
resSig <- resSig[order(resSig$padj,decreasing=FALSE),]
write.csv(resSig,"reports/subset1_Normox_vs_Hypoxia.csv")
write.csv(fpm(dds),"reports/subset1_FPM.csv")
write.csv(fpkm(dds),"reports/subset1_FPKM.csv")

dds$id <- factor(paste0(dds$genotype,".",dds$condition))
dds.Coll <- collapseReplicates(dds,dds$id)
colData(dds.Coll)
vsd.Coll <- vst(dds.Coll, blind=FALSE)

dds.Coll <- DESeq(dds.Coll)
res.Coll <- results(dds.Coll)

mcols(dds.Coll)$basepairs = gene_lengths[rownames(mcols(dds.Coll, use.names=TRUE))]

resSig.Coll <- subset(res.Coll,res.Coll$padj < 0.01)
resSig.Coll <- resSig.Coll[order(resSig.Coll$padj,decreasing=FALSE),]

write.csv(resSig.Coll,"reports/subset1_Collapsed_Normox_vs_Hypoxia.csv")
write.csv(fpm(dds.Coll),"reports/subset1_Collapsed_FPM.csv")
write.csv(fpkm(dds.Coll),"reports/subset1_Collapsed_FPKM.csv")

resBest.Coll <- subset(resSig.Coll,abs(resSig.Coll$log2FoldChange) > 2)

mat.Coll <- dds.Coll[ rownames(resBest.Coll), ]
colconditions.Coll = as.data.frame(colData(dds.Coll)[,c("genotype","condition")])

#mat.Coll.reorder <- assay(mat.Coll)[,c("AF293.Normoxia","hrmA_REV.Normoxia",
#                                     "AF293.Hypoxia","hrmA_REV.Hypoxia")]
pheatmap(assay(mat.Coll),method="complete",
         main = "Collapsed Reps, p-value < 0.01 and log_fold_change > 2", 
         show_colnames=T, show_rownames = F,
         annotation_legend = TRUE, legend=T, cluster_rows=TRUE, 
         cluster_cols = F, cexRow=0.3,
         annotation_col=colconditions.Coll)

