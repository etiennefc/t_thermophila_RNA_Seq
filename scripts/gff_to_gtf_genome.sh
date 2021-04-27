#!/usr/bin/bash

gff_path=data/references/genome.gff3
gtf_output_path=data/references/genome.gtf

#Transform GFF to GTF
gffread -F --keep-genes $gff_path > ${gtf_output_path}_temp

#Change mRNA feature to transcript, add '"' at the end of each line and remove CDS lines in the gtf
sed -i 's/mRNA/transcript/g; s/$/"/g' ${gtf_output_path}_temp
grep -v CDS ${gtf_output_path}_temp > $gtf_output_path
rm ${gtf_output_path}_temp

#Modify last column of the GTF file to have the correct features used by CoCo
sed -i 's/;/"; /g; s/=/ "/g; s/Name/gene_name/g; s/""/"/g' $gtf_output_path #add '"' on right and left of feature, remove '=', remove two consecutive '"'
sed -i 's/Note "/gene_version "1"; gene_biotype "protein_coding"; gene_source "AUGUSTUS"; note "/g; s/geneID/gene_id/g' $gtf_output_path #add gene version, biotype and source; change geneID for gene_id
sed -i 's/ID/transcript_id/g' $gtf_output_path #change ID for transcript_id
sed -i -E 's/transcript_id\s"g[0-9]*";\s//g' $gtf_output_path #remove transcript_id from gene lines
sed -i 's/\("g[0-9]*.t[0-9]*";\)/\1 transcript_name \1 transcript_version "1"; transcript_source "AUGUSTUS"; transcript_support_level "NA"; tag "basic"; transcript_biotype "protein_coding";/g' $gtf_output_path #locate transcript_id, copy it in transcript_name and add transcript version, source, suport level, biotype and tag
sed -i 's/Parent\s\("g[0-9]*";*\)/gene_id \1/g' $gtf_output_path #change Parent for gene_id in transcript lines

gawk -i inplace 'BEGIN{exon_nb=0;id=""}{if (match($0, /(Parent\s"g[0-9]*.t[0-9]*")/, pattern)) {if (pattern[1]==id){exon_nb+=1} else {exon_nb=1} print $0 "; exon_number \""exon_nb"\";"; id=pattern[1]} else print $0}' $gtf_output_path #add exon_number to end of exon lines
gawk -i inplace 'BEGIN{all_features="";id=""}{if (match($0, /(transcript_id\s"g[0-9]*.t[0-9]*"[^\n]*)/, features)) {split(features[1], pattern, " "); all_features=features[1]; id=pattern[2]}; if ($0 ~ "Parent "id){print $0,all_features} else print $0}' $gtf_output_path #add gene and transcript features to exon lines
sed -i -E '/^#/d; s/$/;/g; s/Parent\s("g[0-9]*).t[0-9]*"/exon_id \1"; exon_version "1"/g' $gtf_output_path #remove lines starting with "#", add ';' at the end of each line, remove Parent feature in exon lines and replace it with exon_id (which is the gene_id) and exon_version
