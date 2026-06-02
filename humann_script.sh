#!/bin/bash -l
#SBATCH --job-name=run_humann
#SBATCH --account=project_2009650
#SBATCH --partition=small
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt

module load humann

for file in ./processed_fastq/concat_*.fastq
do
  humann --taxonomic-profile bracken_to_metaphlan/${file:25:3}_metaphlan.txt --threads=$SLURM_CPUS_PER_TASK -i $file --nucleotide-database ../HUMANN/HUMANN_NUC/ --bypass-translated-search  --metaphlan-options "--offline --bowtie2db $MPA" -o hmp_subset
done
