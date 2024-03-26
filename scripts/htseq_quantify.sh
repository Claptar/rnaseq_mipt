#!/bin/bash

mkdir -p countsHtseq

for infile in sortedbam/*.bam
do
  sample=$(basename -a -s _sorted.bam ${infile})
   htseq-count --format bam \
   --stranded no \
   -i gene_id \
   --counts_outpu countsHtseq/${sample}_counts.tsv \
   $infile \
   references/genome.gtf
done