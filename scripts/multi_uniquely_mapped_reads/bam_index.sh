#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=32000M
#SBATCH --array=[0-8]

module load StdEnv/2020
module load samtools/1.12

namelist=(\
'KO_1' \
'KO_2' \
'KO_3' \
'MLP1_IP_1' \
'MLP1_IP_2' \
'MLP1_IP_3' \
'WT_1' \
'WT_2' \
'WT_3' \
)

name=${namelist[$SLURM_ARRAY_TASK_ID]}
bam_path=$SCRATCH/t_thermophila_rna_seq/results/star

samtools index $bam_path/$name/*.bam $bam_path/$name/Aligned.sortedByCoord.out.bam.bai






