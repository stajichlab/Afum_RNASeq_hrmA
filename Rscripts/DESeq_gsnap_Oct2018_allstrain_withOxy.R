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
#samples <- read.csv("samples_plotsubset1.csv",header=TRUE)
samples <- read.csv("samples_nobadrep.csv",header=TRUE)
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
#rld.oxy <- rlog(dds.oxy, blind=FALSE)
head(assay(vsd.oxy), 3)
#head(assay(rld.oxy), 3)

# test
dds.oxy <- DESeq(dds.oxy)
res.oxy <- results(dds.oxy)

FPKM <- fpkm(dds.oxy)

min_profile <- data.frame(apply(FPKM, 1, min) )
rownames(min_profile) <- rownames(FPKM)
colnames(min_profile) <- c("Value")
head(min_profile)
filter_genes <- subset(min_profile,min_profile$Value >= 5)
filter_res.oxy <- subset(res.oxy,rownames(res.oxy) %in% rownames(filter_genes))
nrow(filter_res.oxy)

select <- order(filter_res.oxy$padj,decreasing=TRUE)[1:75]

df2 <- as.data.frame(colData(dds.oxy)[,c("condition","genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Condition","Genotype")
pdf("plots/heatmaps_Oxygenalso_allstrains.pdf")

pheatmap(assay(vsd.oxy)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Oxygen VSD Top Expression")

# end test

res.oxy <- results(dds.oxy, contrast=c("condition","Normoxia","Hypoxia"))

resLFC.oxy <- lfcShrink(dds.oxy, coef=2, type="apeglm")
summary(resLFC.oxy)

# Get diff expressed
resLFC.oxy <- subset(resLFC.oxy,rownames(resLFC.oxy) %in% rownames(filter_genes))

resSig <- subset(resLFC.oxy, resLFC.oxy$padj < 0.01 & 
                   abs(resLFC.oxy$log2FoldChange) >= 2)
resSig <- resSig[order(resSig$padj,decreasing=FALSE),]
write.csv(resSig,"reports/subset1_condition_Normoxia_vs_Hypoxia.csv")
write.csv(fpm(dds.oxy),"reports/subset1_FPM.csv")
write.csv(fpkm(dds.oxy),"reports/subset1_FPKM.csv")

pheatmap(assay(vsd.oxy[rownames(resSig),]), cluster_rows=TRUE,
					     show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,  
         cluster_cols=FALSE, annotation_col=df2,main="VSD p<0.01 and 2-fold Exp Fold change Oxy")

pheatmap(assay(vsd.oxy[rownames(resSig),]), cluster_rows=TRUE,
         show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,  
         cluster_cols=TRUE, annotation_col=df2,main="VSD p<0.01 and 2-fold Exp Fold change Oxy")

pheatmap(assay(vsd.oxy[rownames(resSig)[1:75],]), cluster_rows=TRUE,
         show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,  
         cluster_cols=FALSE, annotation_col=df2,main="Top 50 VSD p<0.01 and 2-fold Exp Fold change Oxy")

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
#rld.geno <- rlog(dds.geno, blind=FALSE)
head(assay(vsd.geno), 3)
#head(assay(rld.geno), 3)

# test
dds.geno <- DESeq(dds.geno)
res.geno <- results(dds.geno)

filter_res.geno <- subset(res.geno,rownames(res.geno) %in% rownames(filter_genes))

select <- order(filter_res.geno$padj,
                decreasing=TRUE)

df2 <- as.data.frame(colData(dds.geno)[,c("condition","genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Condition","Genotype")
pheatmap(assay(vsd.geno)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Genotype VSD All Expression Diff")

pheatmap(assay(vsd.geno)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=TRUE, annotation_col=df2,main="Genotype VSD All Expression Diff")
