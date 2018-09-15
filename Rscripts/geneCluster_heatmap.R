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
library(EDASeq)

# how much wider to make gene cluster range when drawing
EXPAND_GENE_RANGE = 5

my_pal2 = mypal2 <- colorRampPalette(brewer.pal(6, "YlOrRd"))

countdata <- read.table("reports/Hypoxia.Af293.gsnap_reads.nostrand.tab", header=TRUE, row.names=1)
colnames(countdata) <- gsub("\\.gsnap_Afum_Af293\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("aln\\.", "", colnames(countdata))
countdata <- countdata[ ,6:ncol(countdata)]
countdata <- as.matrix(countdata)

samples <- read.csv("samples.csv",header=TRUE)
exprnames <- do.call(paste,c(samples[c("Strain","Condition","Replicate")],sep="_"))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)

# check that experimental columns match in order
#all(exprnames %in% colnames(countdata))
#all(exprnames == colnames(countdata))

# reorder the columns anyways... in case data change along the way
countdata <- countdata[,exprnames]
all(exprnames == colnames(countdata))

pdf("plots/hrmA_Cluster_heatmap.pdf")
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

dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

dds <- DESeq(dds)
res <- results(dds)

#resLFC <- lfcShrink(dds, coef="condition_Normoxia_vs_Hypoxia", type="apeglm")

# plot cluster heatmap/range
seqinfo = 
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
ClusterSelect <- match(ClusterGeneNames,rownames(rld))
pheatmap(assay(rld)[ClusterSelect,],method="complete",main = "hrmA cluster", 
         show_rownames = T,
         annotation_legend = FALSE, legend=T, cluster_rows=FALSE, cexRow=0.55 )


pdf("plots/TopDiffExpFoldChange.pdf")
resTop <- subset(res,abs(res$log2FoldChange) > 4 & res$padj < 0.01)
topVar <- head(resTop[order(resTop$log2FoldChange,decreasing=FALSE),],60)
mat  <- assay(rld)[ rownames(topVar), ]

# reorder the columns
pheatmap(mat,method="complete",main = "TopVar RLD order by fold change", show_rownames = T, show_colnames=T,
         annotation_legend = FALSE, legend=T, cluster_rows=FALSE, cexRow=0.3,
         fontsize_row = 7,fontsize_col = 7)
pheatmap(mat,method="complete",main = "TopVar RLD cluster rows", show_rownames = T, show_colnames=T,
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.3,
         fontsize_row = 7,fontsize_col = 7)

dev.off()

