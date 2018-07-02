#!/bin/bash
#SBATCH --ntasks 1 --nodes 1 --mem 16gb  -p short --out logs/QoRTs.%a.log

GENOME=Afum_Af293
GENOMEFA=genome/AfumigatusAf293_Genome.fasta.gz
CHROMSIZES=genome/AfumigatusAf293_Genome.chromsizes
GFF=genome/Afumigatus_Af293.gtf
FQDIR=data/180109_DA003
INDIR=aln
OUTDIR=results/gsnap_QoRTs
SAMPLEFILE=samples.csv

mkdir -p $OUTDIR

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "cannot run without a number provided either cmdline or --array in sbatch"
 exit
fi

MEM=16
IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read EXP SAMPLE CONDITION REP
do 
 INFILE=$INDIR/${SAMPLE}_${CONDITION}.r${REP}.gsnap_${GENOME}.bam
 OUTDIRRUN=$OUTDIR/${SAMPLE}_${CONDITION}.r${REP} 
 READS=$(ls $FQDIR/${EXP}_*.fastq.gz | perl -p -e 's/\s+/,/g' | perl -p -e 's/,$//')
 echo "$EXP $SAMPLE $CONDITION $REP $INFILE"
  echo "#!/bin/bash" > job_qorts_$EXP.sh 
  echo "module load java/8" >> job_qorts_$EXP.sh
  echo "module load QoRTs" >> job_qorts_$EXP.sh
 if [ ! -d $OUTDIRRUN ]; then
  echo "java -Xmx${MEM}G -jar \$QORTS QC --singleEnded --minMAPQ 20 --stranded  --genomeFA $GENOMEFA \
    --rawfastq $READS \
    $INFILE $GFF $OUTDIRRUN" >> job_qorts_$EXP.sh
 fi
 if [ ! -f $OUTDIRRUN/QC.wiggle.fwd.wig.gz ]; then
  echo "java -Xmx${MEM}G -jar \$QORTS bamToWiggle \
   --singleEnded --stranded \
   --negativeReverseStrand \
   --includeTrackDefLine \
   --minMAPQ 20 \
   $INFILE ${SAMPLE}_${CONDITION}.r${REP} \
   $CHROMSIZES \
   $OUTDIRRUN/QC.wiggle" >> job_qorts_${EXP}.sh
 fi
 bash job_qorts_$EXP.sh
 unlink job_qorts_$EXP.sh
done
