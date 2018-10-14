# compute rna-splicing event differentiation?

decoder.file <- "decoder.sample_Id.txt"
decoder <- read.table(decoder.file, header=T, stringsAsFactors=F)

samples = decoder$sample.ID
Strains = decoder$Strain
Treatment = decoder$Condition

samples
Strains
Treatment

directory <- "results/gsnap_QoRTs_merged/"

countFiles <- paste0(directory,samples,
                     "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz")
countFiles                          
library(JunctionSeq)

design <- data.frame(condition = factor(Treatment),
                     genotype  = factor(Strains) );


#gff.file <- "results/withNovel.forJunctionSeq.gff.gz"
gff.file <- "genome/JunctionSeq.flat.gff.gz"

jscs <- readJunctionSeqCounts(countfiles = countFiles,
                               samplenames = samples,
                               design = design,
                               flat.gff.file = gff.file,
                               nCores = 1,
                               analysis.type = "junctionsAndExons"
);

