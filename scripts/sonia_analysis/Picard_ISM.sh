#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=128000M
#SBATCH --array=[0-9]

nixpkgs/16.09
module load picard/2.18.9

namelist=('decoy' \
'KO_1' \
'KO_2' \
'KO_3' \
'WT_1' \
'WT_2' \
'WT_3' \
'MLP1_IP_1' \
'MLP1_IP_2' \
'MLP1_IP_3' \
)

name=${namelist[$SLURM_ARRAY_TASK_ID]}
project_path=$SCRATCH/t_thermophila_rna_seq/results/

mkdir -p $project_path/Picard/$name

BAM_FILE=$project_path/star/$name/Aligned.sortedByCoord.out.bam

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
      I=$BAM_FILE \
      O=$project_path/Picard/$name/Picard_insert_size_metrics.txt \
      H=$project_path/Picard/$name/Picard_size_histogram.pdf \

