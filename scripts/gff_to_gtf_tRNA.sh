#!/usr/bin/bash

gff_path=data/references/tRNA.gff3
gtf_output_path=data/references/tRNA.gtf

#Transform GFF to GTF
gffread -F --keep-genes $gff_path > $gtf_output_path

#Modify last column of the GTF file to have the correct features used by CoCo
sed -i 's/Note=%/note "/g; s/%/_/g; s/$/"/g; s/;/"; /g; s/=/ "/g; s/ID/gene_id/g; s/Name/gene_name/g' $gtf_output_path

#Duplicate the "tRNA" (gene) lines not the exon lines and substitute the "tRNA" feature for a "transcript" feature in the newly added lines
awk -i inplace 'NR % 2 == 0 {for(i=0;i<1;++i)print}{sub("tRNA", "transcript"); print $0}' $gtf_output_path

#Change remaining tRNA features for "gene" features
sed -i 's/\ttRNA/\tgene/g' $gtf_output_path

#Add gene version, gene biotype and gene source to gene and transcript lines
sed -i 's/note/gene_version "1"; gene_biotype "tRNA"; gene_source "MANUAL_ENTRY"; note/g' $gtf_output_path

#Add transcript id, name, version, source, support level tag and biotype to transcript lines
sed -i -E 's/(transcript.*gene_id ".*"; gene_name)/\1\1/g' $gtf_output_path #duplicate characters in transcript lines contained between transcript and gene_name
sed -i -E 's/gene_id "(.*)"; gene_nametranscript.*gene_id/transcript_id "\1.t1"; transcript_name "\1.t1"; transcript_version "1"; transcript_source "MANUAL_ENTRY"; transcript_support_level "NA"; tag "basic"; transcript_biotype "tRNA"; gene_id/g' $gtf_output_path #add transcript info in transcript lines

#Add gene and transcript info to exon lines
gawk -i inplace 'BEGIN{all_features=""}{if (match($0, /(transcript_id\s"[0-9]*.t[0-9]"[^\n]*)/, features)) {split(features[1], pattern, " "); all_features=features[1]}; if ($0 ~ "Parent "){print $0,all_features} else print $0}' $gtf_output_path

#Remove the first four header lines and add exon info by substituting the Parent feature for exon features
sed -i -E '/^#/d; s/$/;/g; s/Parent "(.*)"\stranscript_id/exon_id "\1"; exon_version "1"; exon_number "1"; transcript_id/g' $gtf_output_path #remove the first four header lines starting with '#', add ';' at the end of each line and add exon_id, exon_version and exon_number to exon lines
