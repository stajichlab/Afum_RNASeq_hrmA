# compute rna-splicing event differentiation?

#library(JunctionSeq)
library(QoRTs)

decoder.file <- "decoder.txt"
directory <- "results/gsnap_QoRTs/"
decoder.data <- read.table(decoder.file, header=T, stringsAsFactors=F)

res <- read.qc.results.data(directory, decoder = decoder.data,
                            calc.DESeq2 = TRUE, calc.edgeR = TRUE)


sizeFactors <- res@calc.data[["norm.factors.bySample"]];
sizeFactors.GEO <- data.frame(sample.ID = sizeFactors$sample.ID,
size.factor = sizeFactors$Norm_Geo);

write.table(sizeFactors.GEO,file = "genome/sizeFactors.GEO.txt",
            row.names=F,col.names=TRUE,sep="\t",quote=F);

makeMultiPlot.all(res, outfile.dir = "./reports/gsnap_QoRTs/",plot.device.name = "pdf");

countFiles <- system.file(paste0(directory,decoder.data$sample.ID,
                                 "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
                          )

gff.file <- "results/gsnap_QoRTs_merged/withNovel.forJunctionSeq.gff.gz"
jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder.data$sample.ID,
                               condition=factor(decoder.data$group.ID),
                               flat.gff.file = gff.file,
                               nCores = 8,
                               analysis.type = "junctionsAndExons"
);

