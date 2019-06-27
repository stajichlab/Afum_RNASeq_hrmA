RNA Seq analysis of eefA A.fumigatus data
Data analysis supports Kowalski et al 2019

Pipeline steps that were used in final manuscript reported data are

pipeline/
  - 00_download_GO.sh - download GOA file
  - 00_download.sh - download genomic and GTF data from FungiDB
  - 00_build_gmap_index.sh - build index files for GMAP/GSNAP run
  - 01_gsnap.sh - run gsnap aligner, run with arrayjobs `sbatch --array=1-30`
  - 02_gsnap_makebam.sh - convert sam to bam, run with arrayjobs `sbatch --array=1-30`
  - 03_gsnap_subread_counts.sh run with arrayjobs `sbatch --array=1-30`

R analysis for processing counts and generating plots were done with 
Rscript Rscript/DESeq_gsnap_Oct2018_allstrains.R
and
Rscript Rscript/DESeq_gsnap_Oct2018_allstrain_withOxy.R

