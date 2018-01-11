#!/bin/bash
#SBATCH --ntasks 1 --nodes 1 --time 1:00:00 -p short -J downloadGO --out logs/GO.download.log

mkdir -p GO
cd GO
RELEASE=35
if [ ! -f 22118.N_fumigata_ATCC_MYA-4609.goa ]; then
 curl -O ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/22118.N_fumigata_ATCC_MYA-4609.goa
fi
if [ ! -f FungiDB-${RELEASE}_Afumigatus_Af293_InterproDomains.txt ]; then
 curl -O http://fungidb.org/common/downloads/release-${RELEASE}/Afumigatus_Af293/txt/FungiDB-${RELEASE}_AfumigatusAf293_InterproDomains.txt
fi
cd ..
perl scripts/goaAfum2simpleGO.pl GO/22118.N_fumigata_ATCC_MYA-4609.goa > GO/Afum.GO.tab
