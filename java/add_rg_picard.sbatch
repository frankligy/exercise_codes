#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=200Gb
#SBATCH --job-name="picard"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err

module load jdk/17u028
java -version

cd bwa_mem2_dna_bam
java -jar ../picard.jar AddOrReplaceReadGroups \
    I="052_7-N_align.bam" \
    O="052_7-N_align_with_RG.bam" \
    SORT_ORDER=coordinate \
    RGID=foo \
    RGLB=bar \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM="052_7-N_align" \
    CREATE_INDEX=True
