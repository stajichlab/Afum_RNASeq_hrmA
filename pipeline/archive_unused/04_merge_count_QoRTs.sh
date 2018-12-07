#!/bin/bash
#SBATCH --ntasks 2 --nodes 1 --mem 16gb  -p short --out logs/QoRTs_mergeCounts.log

GENOME=Afum_Af293
GENOMEFA=genome/AfumigatusAf293_Genome.fasta.gz
INDIR=results/gsnap_QoRTs/
OUTDIR=results/gsnap_QoRTs_merged/
SAMPLEFILE=decoder.txt
SIZEFACTORS=genome/sizeFactors.GEO.txt
JUNCTIONGFF=genome/JunctionSeq.flat.gff.gz
ANNOTGFF=genome/Afumigatus_Af293.gtf.gz

MEM=16
if [ ! -d $OUTDIR ]; then
 mkdir -p $OUTDIR
 java -Xmx${MEM}G -jar $QORTS mergeAllCounts $INDIR $SAMPLEFILE $OUTDIR
fi

for sample in $(tail -n +2 $SAMPLEFILE | cut -f2 | sort | uniq)
do
 if [ ! -f $OUTDIR/$sample/QC.junctionBed.known.bed.gz ]; then
  java -Xmx${MEM}G -jar $QORTS \
   makeJunctionTrack \
   --nonflatgtf \
   --stranded \
   --filenames $OUTDIR/$sample/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz \
   $GFF $OUTDIR/$sample/QC.junctionBed.known.bed.gz
 fi
done

java -Xmx${MEM}G -jar $QORTS \
 mergeNovelSplices \
 --minCount 10 \
 --stranded \
 $OUTDIR $SIZEFACTORS $ANNOTGFF $OUTDIR
