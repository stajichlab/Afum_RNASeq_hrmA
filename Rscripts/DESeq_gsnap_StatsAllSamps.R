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

dds.oxy <- estimateSizeFactors(dds.oxy)
dds.oxy <- estimateDispersions(dds.oxy)

vsd.oxy <- vst(dds.oxy, blind=FALSE)
rld.oxy <- rlog(dds.oxy, blind=FALSE)
head(assay(vsd.oxy), 3)
head(assay(rld.oxy), 3)

# test
dds.oxy <- DESeq(dds.oxy)
res.oxy <- results(dds.oxy)
res.oxy <- results(dds.oxy, contrast=c("condition","Normoxia","Hypoxia"))

resLFC.oxy <- lfcShrink(dds.oxy, coef=2, type="apeglm")
summary(resLFC.oxy)

# Get diff expressed
resSig <- subset(resLFC.oxy, resLFC.oxy$padj < 0.05 & 
                   abs(resLFC.oxy$log2FoldChange) > 2)
resSig <- resSig[order(resSig$padj,decreasing=FALSE),]
write.csv(resSig,"reports/all_condition_Normoxia_vs_Hypoxia.csv")
write.csv(fpm(dds.oxy),"reports/all_FPM.csv")
write.csv(fpkm(dds.oxy),"reports/all_FPKM.csv")


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

dds.geno <- DESeq(dds.geno)
res.geno <- results(dds.geno)
res.geno <- results(dds.geno, contrast=c("genotype","AF293","hrmA_REV"))

resLFC.geno <- lfcShrink(dds.geno, coef=2, type="apeglm")
summary(resLFC.geno)

# Get diff expressed
resSig <- subset(resLFC.geno, resLFC.geno$padj < 0.05 & 
                   abs(resLFC.geno$log2FoldChange) > 2)
resSig <- resSig[order(resSig$padj,decreasing=FALSE),]
write.csv(resSig,"reports/all_genotype_AF293_vs_hrmA_REV.csv")

# == now generate collapsed
dds.oxy$id <- factor(paste0(dds.oxy$genotype,'.',dds.oxy$condition))

dds.oxy.Coll <- collapseReplicates(dds.oxy,dds.oxy$id)
colData(dds.oxy.Coll)

dds.oxy.Coll <- DESeq(dds.oxy.Coll)
res.oxy.Coll <- results(dds.oxy.Coll)

mcols(dds.oxy.Coll)$basepairs = gene_lengths[rownames(mcols(dds.oxy.Coll, use.names=TRUE))]

resSig.oxy.Coll <- subset(res.oxy.Coll,res.oxy$padj < 0.05)
resSig.oxy.Coll <- resSig.oxy.Coll[order(resSig.oxy.Coll$pvalue,decreasing=FALSE),]

write.csv(resSig.oxy.Coll,"reports/all_Collapsed_Normox_vs_Hypoxia.csv")
write.csv(fpm(dds.oxy.Coll),"reports/all_Collapsed_FPM.csv")
write.csv(fpkm(dds.oxy.Coll),"reports/all_Collapsed_FPKM.csv")

# == dds geno
dds.geno$id <- factor(paste0(dds.geno$genotype,'.',dds.geno$condition))

dds.geno.Coll <- collapseReplicates(dds.geno,dds.geno$id)
colData(dds.geno.Coll)

dds.geno.Coll <- DESeq(dds.geno.Coll)
res.geno.Coll <- results(dds.geno.Coll)

mcols(dds.geno.Coll)$basepairs = gene_lengths[rownames(mcols(dds.geno.Coll, use.names=TRUE))]

resSig.geno.Coll <- subset(res.geno.Coll,res.geno$padj < 0.05)
resSig.geno.Coll <- resSig.geno.Coll[order(resSig.geno.Coll$padj,decreasing=FALSE),]
write.csv(resSig.geno.Coll,"reports/subset1_Collapsed_AF293_vs_hrmA_REV.csv")

