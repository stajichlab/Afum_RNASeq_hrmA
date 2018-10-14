library(GenomicRanges)
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
# strain ordering should be
# WT, delta, oe strain - OE-hmrA, rev, evol

# how much wider to make gene cluster range when drawing
EXPAND_GENE_RANGE = 5

my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

countdata <- read.table("reports/Hypoxia.Af293.gsnap_reads.nostrand.tab", header=TRUE, row.names=1)
colnames(countdata) <- gsub("\\.gsnap_Afum_Af293\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("aln\\.", "", colnames(countdata))

colnames(countdata) <- sub("eefA","hrmA",colnames(countdata),perl=TRUE)
countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)

samples <- read.csv("samples.csv",header=TRUE)
samples$Strain <-  sub("eefA","hrmA",samples$Strain,perl=TRUE)
exprnames <- do.call(paste,c(samples[c("Strain","Condition","Replicate")],sep="_"))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)
exprnames <- sub("eefA","hrmA",exprnames,perl=TRUE)

# check that experimental columns match in order
#all(exprnames %in% colnames(countdata))
#all(exprnames == colnames(countdata))

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
                              design    = ~ genotype + condition )

dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)

dds <- DESeq(dds)
res <- results(dds)

# plot cluster heatmap/range
txdb <- makeTxDbFromGFF("genome/FungiDB-35_AfumigatusAf293.gff",
                        dataSource="FungiDB-35",
                        organism="Aspergillus fumigatus")
saveDb(txdb,"AfumAf293.FungiDB35.db")

leftcluster = list(tx_name="Afu5g14865-T")
rightcluster = list(tx_name="Afu5g14920-T")

left <- transcripts(txdb,filter=leftcluster,use.names=TRUE)
right <- transcripts(txdb,filter=rightcluster,use.names=TRUE)

clusterrange <- GRanges( seqnames(left),
                         IRanges(start(range(left)),
                                 end(range(right))))


clusterTx <- transcriptsByOverlaps(txdb,clusterrange)
clusterTx <- sort(clusterTx, decreasing=FALSE, ignore.strand=TRUE)

clusterTx

ChromTxs = transcripts(txdb,filter=list(tx_chrom = seqnames(left)))
ChromTxs = sort(ChromTxs, decreasing=FALSE, ignore.strand=TRUE)

leftindex = match(left,ChromTxs)
rightindex = match(right,ChromTxs)

left_extend = ChromTxs[leftindex - EXPAND_GENE_RANGE]
right_extend = ChromTxs[rightindex + EXPAND_GENE_RANGE]

clusterExtRange <- GRanges( seqnames(left),
                        IRanges(start(range(left_extend)),
                                end(range(right_extend))))

clusterTxExtended <- transcriptsByOverlaps(txdb,clusterExtRange)
clusterTxExtended <- sort(clusterTxExtended, decreasing=FALSE, ignore.strand=TRUE)

clusterTxExtended

ClusterGeneNames <- sub('-T','',clusterTxExtended$tx_name)
ClusterSelect <- match(ClusterGeneNames,rownames(vsd))

pdf("plots/hrmA_Cluster_heatmap.pdf")

colconditions = as.data.frame(colData(vsd)[,c("genotype","condition")])
pheatmap(assay(vsd)[ClusterSelect,],method="complete",main = "hrmA cluster", 
         show_rownames = T,show_colnames = F,
         annotation_legend = TRUE, legend=T, annotation_col=colconditions,
         cluster_cols=FALSE,cluster_rows=FALSE, cexRow=0.55 )

pheatmap(assay(vsd)[ClusterSelect,],method="complete",main = "hrmA cluster", 
         show_rownames = T,show_colnames = T,
         annotation_legend = TRUE, legend=T, annotation_col=colconditions,
         cluster_cols=FALSE,cluster_rows=FALSE, cexRow=0.55 )


#pdf("plots/TopDiffExpFoldChange.pdf")
#resLFC <- lfcShrink(dds, coef="condition_Normoxia_vs_Hypoxia", type="apeglm")
#abs(res$log2FoldChange) > 4 & 
resBest <- subset(res,res$padj < 0.01)
mat  <- assay(vsd)[ rownames(resBest), ]

colconditions = as.data.frame(colData(vsd)[,c("genotype","condition")])
pheatmap(mat,method="complete",main = "HeatMap of Expression (p-value < 0.01)", show_colnames=F, show_rownames = F,
         annotation_legend = TRUE, legend=T, cluster_rows=TRUE, cluster_cols = F, cexRow=0.3,
         annotation_col=colconditions)

topVar <- head(resBest[order(resBest$log2FoldChange,decreasing=FALSE),],100)
mat  <- assay(vsd)[ rownames(topVar), ]
pheatmap(mat,method="complete",main = "Top 100 based on p-value", show_colnames=F, show_rownames = T,
         annotation_legend = TRUE, legend=T, cluster_rows=TRUE, cluster_cols = F, cexRow=0.3,
         fontsize_row = 6,fontsize_col = 7,
         annotation_col=colconditions)

dds$id <- factor(paste0(dds$genotype,".",dds$condition))
ddsColl <- collapseReplicates(dds,dds$id)

colData(ddsColl)
#ddsColl <- estimateSizeFactors(ddsColl)
#ddsColl <- estimateDispersions(ddsColl)
vsdColl <- vst(ddsColl, blind=FALSE)

vsdColl.reorder <- assay(vsdColl)[,c("AF293.Normoxia","AF293.Hypoxia",
                                     "Delta_hrmA_AF293.Normoxia","Delta_hrmA_AF293.Hypoxia",
                                     "hrmA_OE.Normoxia", "hrmA_OE.Hypoxia", 
                                     "hrmA_REV.Normoxia", "hrmA_REV.Hypoxia",
                                     "EVOL.Normoxia","EVOL.Hypoxia")]
#dds.Coll <- DESeq(ddsColl)
#res.Coll <- results(dds.Coll)

#pdf("plots/Reps_Collapsed.pdf")
#resLFC.Coll <- lfcShrink(ddsColl, coef="condition_Normoxia_vs_Hypoxia", type="apeglm")
resBest.Coll <- subset(vsdColl.reorder,res$padj < 0.01)
mat.Coll <- vsdColl.reorder[ rownames(resBest.Coll), ]
colconditions.Coll = as.data.frame(colData(vsdColl)[,c("genotype","condition")])

pheatmap(mat.Coll,method="complete",main = "Collapsed Reps, p-value < 0.01", show_colnames=T, show_rownames = F,
         annotation_legend = TRUE, legend=T, cluster_rows=TRUE, cluster_cols = F, cexRow=0.3,
         annotation_col=colconditions.Coll)

ClusterSelect.Coll <- match(ClusterGeneNames,rownames(vsdColl))

pheatmap(vsdColl.reorder[ClusterSelect.Coll,],method="complete",main = "hrmA cluster", 
         show_rownames = T,
         legend=TRUE, 
         cluster_cols=FALSE,cluster_rows=FALSE, cexRow=0.55, 
         annotation_col = colconditions.Coll )


clustermat.Coll  <- vsdColl.reorder[ ClusterSelect, ]
clustermat.Coll.RowMeans  <- clustermat.Coll - rowMeans(clustermat.Coll)
pheatmap(clustermat.Coll.RowMeans,method="complete",main = "hrmA cluster rowMean normalized", 
         show_rownames = T,
         legend=TRUE, 
         cluster_cols=FALSE,cluster_rows=FALSE, cexRow=0.55, annotation_col = colconditions.Coll )


clustermat.Coll.Af293Norm  <- log(clustermat.Coll / clustermat.Coll[,c("AF293.Normoxia")]) / log(2)
pheatmap(clustermat.Coll.Af293Norm,method="complete",main = "hrmA cluster Af293 normalized", 
         show_rownames = T,
         legend=TRUE, 
         cluster_cols=FALSE,cluster_rows=FALSE, cexRow=0.55,annotation_col = colconditions.Coll  )

# now whole dataset

mat.Coll.RowMeans  <- vsdColl.reorder - rowMeans(vsdColl.reorder)
pheatmap(mat.Coll.RowMeans,method="complete",main = "Pvalue less than 0.01 rowMeans normalized", 
         show_rownames = F,
         legend=TRUE, 
         cluster_cols=FALSE,cluster_rows=TRUE, cexRow=0.55,annotation_col = colconditions.Coll  )

mat.Coll.Af293Norm <- log(vsdColl.reorder / vsdColl.reorder[,c("AF293.Normoxia")]) / log(2)

pheatmap(mat.Coll.Af293Norm,method="complete",main = "Pvalue less than 0.01 Af293Norm normalized", 
         show_rownames = F,
         legend=TRUE, 
         cluster_cols=FALSE,cluster_rows=TRUE, cexRow=0.55,annotation_col = colconditions.Coll  )
