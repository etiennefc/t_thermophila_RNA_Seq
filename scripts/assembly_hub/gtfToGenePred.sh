#!/bin/bash

module load nixpkgs/16.09  
module load gcc/5.4.0
module load kentutils/20180716

# Convert gtf file into GenePred format, then bigGenePred format (required format for Assembly hub display on UCSC genome browser)

gtf_path=$SCRATCH/t_thermophila_rna_seq/data/references/t_thermophila_complete_2021.gtf
index_path=$SCRATCH/t_thermophila_rna_seq/data/references
output_path=$SCRATCH/t_thermophila_rna_seq/data/references/t_thermophila_complete_2021.bb
chr_size=$SCRATCH/t_thermophila_rna_seq/data/star_index/chrNameLength.txt

gtfToGenePred $gtf_path temp.genePred

genePredToBigGenePred temp.genePred temp.txt 

# The following genes 80406, 80407 are not contained within chr_size file and they are not expressed, so we remove them
LC_COLLATE=C sort -k1,1 -k2,2n temp.txt | sed '/80406/d; /80407/d; /80234/d; /5S.19-1.t1/d; /5S.12-1.t1/d' > temp2.txt

wget https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as

bedToBigBed -extraIndex=name,geneName -type=bed12+8 -tab -as=bigGenePred.as temp2.txt $chr_size $output_path

# Create search index (so we can search by gene name in the Assembly hub)
cat temp.genePred | awk '{print $1, $12, $1}' > temp_input.txt
ixIxx temp_input.txt $index_path/search_index.ix $index_path/search_index.ixx


rm bigGenePred.as
rm temp*
echo "All done!"


