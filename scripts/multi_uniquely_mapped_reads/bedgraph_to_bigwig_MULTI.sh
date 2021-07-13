#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --mem=32000M
#SBATCH --array=[0-8]

module load nixpkgs/16.09  
module load gcc/5.4.0
module load kentutils/20180716



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
bedgraph_path=$SCRATCH/t_thermophila_rna_seq/results/coco_bedgraph_MULTI/
output_path=$SCRATCH/t_thermophila_rna_seq/results/bedgraph_to_bigwig_MULTI/

mkdir -p $output_path


sort -k1,1 -k2,2n $bedgraph_path/${name}.bedgraph > $output_path/${name}_sorted.bedGraph &&


bedGraphToBigWig $output_path/${name}_sorted.bedGraph \
$SCRATCH/t_thermophila_rna_seq/data/star_index/chrNameLength.txt \
$output_path/${name}.bigwig

echo "Well done!"


