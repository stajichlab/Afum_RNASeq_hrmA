#!/usr/bin/bash
#SBATCH --nodes 1 --mem 1G --time 0:30:0

cd genome
RELEASENUM=35
if [ ! -f FungiDB-${RELEASENUM}_AfumigatusAf293.gff ]; then
 curl -O http://fungidb.org/common/downloads/release-$RELEASENUM/AfumigatusAf293/gff/data/FungiDB-${RELEASENUM}_AfumigatusAf293.gff
fi
if [ ! -f FungiDB-${RELEASENUM}_AfumigatusAf293_Genome.fasta ]; then
 curl -O http://fungidb.org/common/downloads/release-$RELEASENUM/AfumigatusAf293/fasta/data/FungiDB-${RELEASENUM}_AfumigatusAf293_Genome.fasta
 if [ -f AfumigatusAf293_Genome.fasta ]; then
  rm AfumigatusAf293_Genome.fasta
 fi
 ln -s FungiDB-${RELEASENUM}_AfumigatusAf293_Genome.fasta AfumigatusAf293_Genome.fasta
fi
if [ ! -f FungiDB-${RELEASENUM}_AfumigatusAf293_AnnotatedCDSs.fasta  ]; then
 curl -O http://fungidb.org/common/downloads/release-$RELEASENUM/AfumigatusAf293/fasta/data/FungiDB-${RELEASENUM}_AfumigatusAf293_AnnotatedCDSs.fasta
 ln -s FungiDB-${RELEASENUM}_AfumigatusAf293_AnnotatedCDSs.fasta AfumigatusAf293_CDS.fasta
fi
grep -P "\t(CDS|exon)\t" FungiDB-${RELEASENUM}_AfumigatusAf293.gff >  Afumigatus_Af293.genes_${RELEASENUM}.gff3
perl -p -e 's/;protein_source_id\S+//; s/ID=[^;]+;Parent=((\S+)-T\S*)(;protein\S+)?/gene_id \"$2\"; transcript_id \"$1\";/; s/ID=[^;]+;//' Afumigatus_Af293.genes_${RELEASENUM}.gff3 > Afumigatus_Af293.v${RELEASENUM}.gtf
ln -s Afumigatus_Af293.v${RELEASENUM}.gtf Afumigatus_Af293.gtf

module load hisat2
hisat2_extract_splice_sites.py Afumigatus_Af293.gtf  > Afumigatus_Af293.v${RELEASENUM}.splicesite.bed
hisat2_extract_exons.py Afumigatus_Af293.gtf  > Afumigatus_Af293.v${RELEASENUM}.exons.bed
ln -s Afumigatus_Af293.v${RELEASENUM}.splicesite.bed Afumigatus_Af293.splicesite.bed
ln -s Afumigatus_Af293.v${RELEASENUM}.exons.bed Afumigatus_Af293.exons.bed

module load kallisto
kallisto index -i Afumigatus_Af293.kallisto.idx AfumigatusAf293_CDS.fasta
