#!/bin/bash

mkdir -p sortedbam

for infile in hisatOutput/*
do
  sample=$(basename -a ${infile})
  samtools view -@4 -b hisatOutput/${sample}/alignment.sam > hisatOutput/${sample}/alignment.bam
  samtools sort -@4 hisatOutput/${sample}/alignment.bam -o sortedbam/${sample}_sorted.bam
  samtools index -@4 sortedbam/${sample}_sorted.bam
done