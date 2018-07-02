#!/usr/bin/bash
#SBATCH --time 2:00:00 --mem 2G 

module load gmap/2018-02-12
cd genome
gmap_build -D=. -d Afum_Af293 AfumigatusAf293_Genome.fasta
gtf_splicesites < Afumigatus_Af293.gtf > Afumigatus_Af293.splicesites.txt
iit_store -o Afum_Af293/Afum_Af293.maps/splicesites.iit < Afumigatus_Af293.splicesites.txt
