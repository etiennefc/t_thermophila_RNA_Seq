#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=32000M
#SBATCH --array=[0-8]

module load StdEnv/2020
module load samtools/1.12

## Separate uniquely mapped reads and multi-mapped reads into two separate bam files

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

echo "Starting uniquely..."
#Get uniquely mapped read (number of hits (NH) = 1) and keep header lines (@SQ, @HD, @PG, @CO)
samtools view -h $bam_path/$name/*out.bam | grep -Ew 'NH:i:1|@SQ|@HD|@PG|@CO' | samtools view -b -o $bam_path/$name/Aligned.sortedByCoord.out.UNIQUE.bam
echo "uniquely done!"

echo "Starting multi..."
#Get multi-mapped read (number of hits (NH) > 1)
samtools view -h $bam_path/$name/*out.bam | grep -v "NH:i:0" | grep -v -w "NH:i:1" | samtools view -b -o $bam_path/$name/Aligned.sortedByCoord.out.MULTI.bam
echo "Multi done!"






