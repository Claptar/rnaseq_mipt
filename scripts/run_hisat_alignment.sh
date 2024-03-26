#!/usr/bin/env bash

# make directories for quantified samples
mkdir -p hisatOutput

# hisat alignment
for infile in fastq_files/*.fastq.gz
do
  base=$(basename -a -s .fastq.gz ${infile})
  mkdir -p hisatOutput/${base}
  hisat2 -p 8 \
  --phred33 \
  -x references/hisatIndex/musgrcm39Index \
  -U $infile \
  -S hisatOutput/${base}/alignment.sam \
  --summary-file hisatOutput/${base}/summary.txt
done