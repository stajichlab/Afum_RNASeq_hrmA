#!/bin/bash 
#SBATCH --nodes 1 --ntasks 8 --mem 12G -J subreadcount -p short
#SBATCH --time 2:00:00 --out logs/subread_count.log

module load subread

# transcript file was updated to recover missing genes
GFF=genome/Afumigatus_Af293.gtf
OUTDIR=reports
INDIR=aln
GENOME=genome/AfumigatusAf293_Genome.fasta
GENOMENAME=Afum_Af293 # so we can easily switch to a different reference
EXTENSION=gsnap_${GENOMENAME}.bam
mkdir -p $OUTDIR
TEMP=/scratch
SAMPLEFILE=samples.csv

CPUS=$SLURM_CPUS_ON_NODE

IFS=,

OUTFILE=$OUTDIR/Hypoxia.Af293.gsnap_reads.nostrand.tab

if [ ! -f $OUTFILE ]; then
    featureCounts -g gene_id -T $CPUS -G $GENOME -s 0 -a $GFF \
        --tmpDir $TEMP \
	-o $OUTFILE -F GTF $INDIR/*.bam
fi

OUTFILE=$OUTDIR/Hypoxia.Af293.gsnap_reads.MM.nostrand.tab
if [ ! -f $OUTFILE ]; then
    featureCounts -g gene_id -T $CPUS -G $GENOME -s 0 -a $GFF \
        --tmpDir $TEMP -M -primary \
	-o $OUTFILE -F GTF $INDIR/*.bam
fi

OUTFILE=$OUTDIR/Hypoxia.Af293.gsnap_frags.nostrand.tab
if [ ! -f $OUTFILE ]; then
    featureCounts -g gene_id -T $CPUS -G $GENOME -s 2 -a $GFF \
        --tmpDir $TEMP \
	-o $OUTFILE -F GTF -p $INDIR/*.bam
fi

