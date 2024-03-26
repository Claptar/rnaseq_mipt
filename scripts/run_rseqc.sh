#!/bin/bash

mkdir -p RSeQCresults/bamstat
mkdir -p RSeQCresults/cliping_profile/
mkdir -p RSeQCresults/delition_profile/
mkdir -p RSeQCresults/insertion_profile/
mkdir -p RSeQCresults/mismatch_profile/
mkdir -p RSeQCresults/duplication_profile/

for infile in sortedbam/*.bam
do
  sample=$(basename -a -s _sorted.bam ${infile})
  echo $sample
  bam_stat.py  -i $infile > RSeQCresults/bamstat/${sample}.txt
  clipping_profile.py -i $infile -s "SE" -o RSeQCresults/cliping_profile/$sample
  deletion_profile.py -i $infile -o RSeQCresults/delition_profile/$sample -l 65
  insertion_profile.py -i $infile -s "SE" -o RSeQCresults/insertion_profile/$sample
  mismatch_profile.py -i $infile -o RSeQCresults/mismatch_profile/$sample -l 65
  read_duplication.py -i $infile -o RSeQCresults/duplication_profile/$sample
done