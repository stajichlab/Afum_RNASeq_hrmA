#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 8gb --time 2:00:00 -p short --out logs/STAR.index.log

module load STAR
cd genome
CPU=8
GENOME=AfumigatusAf293_Genome.fasta
GENOMENAME=Afum_Af293.STAR
GFF=Afumigatus_Af293.gtf

mkdir -p $GENOMENAME
STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $GENOMENAME --genomeFastaFiles $GENOME --sjdbGTFfile $GFF --sjdbOverhang 75
