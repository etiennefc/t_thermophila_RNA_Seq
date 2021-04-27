#!/usr/bin/bash

#Generate protein-coding genes gtf from gff3
./scripts/gff_to_gtf_genome.sh

#Generate tRNA gtf from gff3
./scripts/gff_to_gtf_tRNA.sh

#Generate 5s rRNA gtf from gff3
./scripts/gff_to_gtf_5s_rRNA.sh

