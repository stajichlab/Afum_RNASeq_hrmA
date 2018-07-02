#!/usr/bin/bash
SAMPLEFILE=samples.csv
IFS=,
tail -n +2 $SAMPLEFILE | while read EXP SAMPLE CONDITION REP
do
 echo "ln -s ${SAMPLE}_${CONDITION}.r${REP} ${EXP}"
done
