#!/usr/bin/bash

gff_path=data/references/5s_rRNA.gff3
gtf_output_path=data/references/5s_rRNA.gtf


#Modify gene_id so that it is unique to each gene (add "-1" at every first occurence of that gene id, then add "-2" at the second occurence, then add "-3" ... so that exact copies of genes located on different chromosomes but with the same gene length have different gene_id than in the initial gff3 file
gawk 'BEGIN{copy_nb=0;id=""}{if (match($0, /(5S.[0-9]*;)/, pattern)) {if (pattern[1]==id){copy_nb+=1} else {copy_nb=1} print $0 "ID="pattern[1]"-"copy_nb";Name="pattern[1]"-"copy_nb; id=pattern[1]}}' $gff_path > ${gff_path}_temp
sed -i -E 's/;-/-/g; s/\tID.*ID/\tID/g' ${gff_path}_temp #remove unwanted characters before transforming gff into gtf

#Transform GFF to GTF
gffread -F --keep-genes ${gff_path}_temp > $gtf_output_path
rm ${gff_path}_temp

#Modify last column of the GTF file to have the correct features used by CoCo
sed -i 's/$/"/g; s/;/"; /g; s/=/ "/g; s/ID/gene_id/g; s/Name/gene_name/g' $gtf_output_path

#Duplicate the "5SRNA" (gene) lines not the exon lines and substitute the "5SRNA" feature for a "transcript" feature in the newly added lines
awk -i inplace 'NR % 2 == 0 {for(i=0;i<1;++i)print}{sub("5SRNA", "transcript"); print $0}' $gtf_output_path

#Change remaining 5SRNA features for "gene" features
sed -i 's/\t5SRNA/\tgene/g' $gtf_output_path

#Add gene version, gene biotype and gene source to gene and transcript lines
sed -i -E 's/(gene_name ".*")/\1; gene_version "1"; gene_biotype "rRNA"; gene_source "MANUAL"/g' $gtf_output_path

#Add transcript id, name, version, source, support level tag and biotype to transcript lines
sed -i -E 's/(transcript.*gene_id ".*"; gene_name)/\1\1/g' $gtf_output_path #duplicate characters in transcript lines contained between transcript and gene_name
sed -i -E 's/gene_id "(.*)"; gene_nametranscript.*gene_id/transcript_id "\1.t1"; transcript_name "\1.t1"; transcript_version "1"; transcript_source "MANUAL"; transcript_support_level "NA"; tag "basic"; transcript_biotype "rRNA"; gene_id/g' $gtf_output_path #add transcript info in transcript lines

#Add gene and transcript info to exon lines
gawk -i inplace 'BEGIN{all_features=""}{if (match($0, /(transcript_id "5S.[0-9]*-[0-9]*.t[0-9]"[^\n]*)/, features)) {split(features[1], pattern, " "); all_features=features[1]}; if ($0 ~ "Parent "){print $0,all_features} else print $0}' $gtf_output_path

#Remove the first four header lines and add exon info by substituting the Parent feature for exon features
sed -i -E '/^#/d; s/$/;/g; s/Parent "(.*)"\stranscript_id/exon_id "\1"; exon_version "1"; exon_number "1"; transcript_id/g' $gtf_output_path #remove the first four header lines starting with '#', add ';' at the end of each line and add exon_id, exon_version and exon_number to exon lines
