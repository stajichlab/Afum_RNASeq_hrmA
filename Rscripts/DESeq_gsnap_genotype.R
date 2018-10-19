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

FPKM_MIN=2

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

# == now test for difference when using genotype
dds.geno <- DESeqDataSetFromMatrix(countData = countdata,
                                   colData   = sampleTable, 
                                   design    = ~genotype)

nrow(dds.geno)
dds.geno <- dds.geno[ rowSums(counts(dds.geno)) > 1, ]
nrow(dds.geno)

mcols(dds.geno)$basepairs = gene_lengths[rownames(mcols(dds.geno, use.names=TRUE))]
FPKM <- fpkm(dds.geno)

min_profile <- data.frame(apply(FPKM, 1, min) )
rownames(min_profile) <- rownames(FPKM)
colnames(min_profile) <- c("Value")
head(min_profile)

FPKM_filter_genes <- subset(min_profile,min_profile$Value >= FPKM_MIN)

write.csv(fpm(dds.geno),"reports/subset1_FPM.csv")
write.csv(fpkm(dds.geno),"reports/subset1_FPKM.csv")


dds.geno <- estimateSizeFactors(dds.geno)
dds.geno <- estimateDispersions(dds.geno)

vsd.geno <- vst(dds.geno, blind=FALSE)
rld.geno <- rlog(dds.geno, blind=FALSE)
head(assay(vsd.geno), 3)
head(assay(rld.geno), 3)

# test
dds.geno <- DESeq(dds.geno)
res.geno <- results(dds.geno)
res.geno <- results(dds.geno, contrast=c("genotype","AF293","hrmA_REV"))

filter_res.geno <- subset(res.geno,rownames(res.geno) %in% rownames(FPKM_filter_genes))

# test
select <- order(filter_res.geno$padj,
                decreasing=TRUE)[1:100]

df2 <- as.data.frame(colData(dds.geno)[,c("condition","genotype")])
rownames(df2) = exprnames
colnames(df2) = c("Condition","Genotype")

# don't show vsd anymore

pdf("plots/genotype_heatmaps_alldata.pdf")
#pheatmap(assay(vsd.geno)[select,], cluster_rows=TRUE, show_rownames=TRUE,
#         fontsize_row = 7,fontsize_col = 7,
#         cluster_cols=FALSE, annotation_col=df2,main="Geno VSD Top Expression Diff")
pheatmap(assay(rld.geno)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Geno Top 100 Expression Diff")

# end test
# plot all of te 
pheatmap(assay(rld.geno), cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Geno All Expression")

resLFC.geno <- lfcShrink(dds.geno, coef=2, type="apeglm")
summary(resLFC.geno)

select <- rownames(subset(resLFC.geno, resLFC.geno$padj < 0.01))

write.csv(resLFC.geno,
          "reports/genotype_AF293_vs_hrmA_REV.no_filter.csv")

pheatmap(assay(rld.geno)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="Geno filtered by P < 0.05")

dev.off()

# Get diff expressed
write.csv(subset(resLFC.geno,resLFC.geno$padj < 0.05),
                 "reports/genotype_AF293_vs_hrmA_REV.pvalue_filter.csv")

resSig.geno <- subset(resLFC.geno, rownames(resLFC.geno) %in% rownames(FPKM_filter_genes)
                    & resLFC.geno$padj < 0.05
                    & abs(resLFC.geno$log2FoldChange) >= 2)

resSig.geno <- resSig.geno[order(resSig.geno$padj,decreasing=FALSE),]

write.csv(resSig.geno,"reports/genotype_AF293_vs_hrmA_REV.pvalue_FPKM_Filter.csv")

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

dev.off()
pdf("plots/genotype_filtered_heatmaps.pdf",width=10)

pheatmap(rld.geno.reorder[rownames(resSig.geno),], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main=paste("FPKM >= ",FPKM_MIN," p<0.01, 2-fold Exp Fold change Geno"))



# == COLLAPSE REPLICATES
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

# just use the geno significant set from the replicate data
resSig.geno.Coll <- subset(res.geno.Coll,rownames(res.geno.Coll) %in% rownames(resSig.geno)) 
resSig.geno.Coll <- resSig.geno.Coll[order(resSig.geno.Coll$padj,decreasing=FALSE),]

write.csv(resSig.geno.Coll,"reports/genotype_Collapsed_AF293_vs_hrmA_REV.pvalue_FPKM_Filter.csv")

mat.geno.Coll <- rld.geno.Coll[ rownames(resSig.geno.Coll), ]
colconditions.geno.Coll = as.data.frame(colData(dds.geno.Coll)[,c("condition","genotype")])

#mat.Coll.reorder <- assay(mat.Coll)[,c("AF293.Normoxia","hrmA_REV.Normoxia",
#                                     "AF293.Hypoxia","hrmA_REV.Hypoxia")]
pheatmap(assay(mat.geno.Coll),method="complete",
         main = paste("Collapsed Reps - Geno model, p-value < 0.01 and log_fold_change >= 2 FPKM >=", FPKM_MIN),
         show_colnames=TRUE, show_rownames = TRUE,
         annotation_legend = TRUE, legend=TRUE, cluster_rows=TRUE, 
         cluster_cols = FALSE, cexRow=0.3,
         annotation_col=colconditions.geno.Coll)

