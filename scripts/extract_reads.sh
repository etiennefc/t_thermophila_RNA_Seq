#!/usr/bin/bash

# The 1st argument given to this script is the fasta file of tRNA sequences
# The 2nd argument is the bam input file
# The 3rd argument is the output file
seq_file=$1
bam=$2
output=$3


# Iterate through each tRNA sequences
cat $seq_file | while read line; do
    if [[ $line =~ ^">" ]]; then
        temp_id=${line: -5}  # get gene_id from each id lines (starting with '>')
				position=${line::-6}
				position=${position:1}  # get position of tRNA (to filter bam for reads overlapping this region only)
		else
				echo ">"$temp_id >> $output
				seq=${line: -50}; seq2=${seq:0:49}  # get the last 50 or 49 nt in tRNA sequence
				echo "Mature" >> $output;
				samtools view $bam $position | grep -cE "(${seq}CCA|${seq2}CCA)" >> $output;  # get reads of mature tRNAs ending with CCA
				echo "Premature" >> $output;
				samtools view $bam $position | grep -cE "${seq}T{1,5}" >> $output  # get reads of premature tRNAs (ending with 1 to 5 Ts)
    fi
done
