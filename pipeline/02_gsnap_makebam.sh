#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 16gb --time 2:00:00 -p short -J gsnap.makebam --out logs/gsnap_makebam.%a.log

module load picard
GENOME=Afum_Af293
MEM=16
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=$SLURM_ARRAY_TASK_ID
INDIR=aln
OUTDIR=aln
SAMPLEFILE=samples.csv
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "cannot run without a number provided either cmdline or --array in sbatch"
 exit
fi

IFS=,

tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read EXP SAMPLE CONDITION REP
do
 echo $EXP $SAMPLE $REP
infile=${INDIR}/${SAMPLE}_${CONDITION}.r${REP}.gsnap_${GENOME}.sam
bam=${INDIR}/$(basename $infile .sam).bam

echo "$infile $bam"
if [ ! -f $bam -a ! -z $bam ]; then
 java -Xmx${MEM}g -jar $PICARD CleanSam I=$infile O=${bam}.clean
 java -Xmx${MEM}g -jar $PICARD AddOrReplaceReadGroups SO=coordinate I=${bam}.clean O=${bam} CREATE_INDEX=true TMP_DIR=/scratch/${USER} RGID=$EXP RGSM=${SAMPLE}_${CONDITION}.r${REP} RGPL=NextSeq RGPU=Dartmouth RGLB=$EXP VALIDATION_STRINGENCY=LENIENT
fi
    
#if [ -f $bam ] &&  [ ! -z $bam ]; then
#  rm $infile; touch $infile
#  touch $bam
# fi
done
