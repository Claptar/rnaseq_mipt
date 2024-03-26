#!/bin/bash

# скачаем референсный геном для hisat2 с ENSEMBL
wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
gzip -d Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
mv Mus_musculus.GRCm39.dna_sm.primary_assembly.fa references/genome.fa

# скачаем геномную аннотацию с ENSEMBL
wget https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz
gzip -d Mus_musculus.GRCm39.111.gtf.gz
mv Mus_musculus.GRCm39.111.gtf references/genome.gtf

# скачаем референсный транскриптом для kallisto с ENSEMBL
wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
gzip -d Mus_musculus.GRCm39.cdna.all.fa.gz
mv Mus_musculus.GRCm39.cdna.all.fa references/transcriptome.fa