#!/bin/bash

# make directories for quantified samples
mkdir -p kallistoOutput

# kallisto quantification step
for infile in fastq_files/*.fastq.gz
do
  sample=$(basename -a -s .fastq.gz ${infile})
  echo ${sample}
  kallisto quant -i references/musGRCm39_kallisto -o kallistoOutput/${sample} --single -l 180 -s 20.0 -t 8 fastq_files/${sample}.fastq.gz &> kallistoOutput/${sample}/stdout.log
done