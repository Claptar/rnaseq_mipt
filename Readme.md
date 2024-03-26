# RNA-seq

- Результат
    
    

# Введение

Целью данного задания является сравнение RNA-seq данных перепрограммированных и неперепрограммированных (контрольных) мышиных эмбриональных фибробластов (MEFs) и нахождение генов, которые наиболее сильно изменяют свою экспрессию в этом процессе. В качестве дополнительного задания предлагается изучить основные. 

# Подготовим рабочее пространство

Я работаю на компухтере с операционной системой Windows. Поэтому перед началом работы были установлены `WSL` и `Ubuntu 22.04.3`. Так же дополнительно была установлена `Anaconda3-2023.09-0-Linux-x86_64`

Переда началом подготовим рабочее пространство, а именно создадим рабочую директорию и все необходимые вложенные файлы.

```bash
# Создадим рабочую дирректорию и переместимся в неё
$ mkdir -p omics_hw/rnaseq
$ cd omics_hw/rnaseq

# Cоздадим необходимые дирректории
$ mkdir fastq_files # дирректория с чистыми .fastq файлы
$ mkdir scripts # дирректория со скриптами с помощью которых мы будем обрабатывать данные
$ mkdir references # дирректория c геномным референсом, аннотацией и индексами для генома
$ mkdir qc # дирректория с отчётами о качестве наших данных
```

# Скачивание файлов

Скачаем данные PRJNA836496 при помощи [SRA Explorer](https://sra-explorer.info/#). Для этого создадим скрипт `load_fastq.sh` в дирректории `scripts`. Добавим `$1` перед названием файла чтобы указать необходимую нам директорию.

```bash
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/009/SRR3414629/SRR3414629.fastq.gz -o $1/SRR3414629_GSM2130680_MEF_OSKM_induced_Day7_Scrambled_shRNA_replicate1_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/000/SRR3414630/SRR3414630.fastq.gz -o $1/SRR3414630_GSM2130681_MEF_OSKM_induced_Day7_Scrambled_shRNA_replicate2_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/001/SRR3414631/SRR3414631.fastq.gz -o $1/SRR3414631_GSM2130682_MEF_OSKM_induced_Day7_Scrambled_shRNA_replicate3_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/006/SRR3414636/SRR3414636.fastq.gz -o $1/SRR3414636_GSM2130687_MEF_Day7_Scrambled_shRNA_replicate2_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/005/SRR3414635/SRR3414635.fastq.gz -o $1/SRR3414635_GSM2130686_MEF_Day7_Scrambled_shRNA_replicate1_Mus_musculus_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/007/SRR3414637/SRR3414637.fastq.gz -o $1/SRR3414637_GSM2130688_MEF_Day7_Scrambled_shRNA_replicate3_Mus_musculus_RNA-Seq.fastq.gz
```

Скачаем `.fastq` файлы для нужных нам образцов

```bash
$ bash scripts/load_fastq.sh fastq_files
```

Переименуем файлы чтобы с ними было проще работать. Для этого сделаем небольшой скрипт на `python` 

```python
import os
import glob

dir_name = 'fastq_files'
for filename in glob.glob(f'{dir_name}/*'):
    sra_code = filename.split('/')[1].split('_')[0]
    new_name = f'{dir_name}/{sra_code}.fastq.gz'
    os.rename(filename, new_name)
    print(sra_code)

```

Запустим скрипт

```bash
$ python3 scripts/rename_fastq.py
```

# QC для сырых образцов

Создадим анакондовское окружение и установим  `fastqc`, `multiqc`

```bash
$ conda create -n preprocess -c conda-forge -c bioconda fastqc multiqc
$ conda activate preprocess
```

Теперь сделаем файл с запуском анализа качества прочтений

```bash
#!/bin/bash
### Run QC
fastqc -o qc -t 8 fastq_files/*.fastq.gz

### Run MultiQC
multiqc qc --outdir qc --filename multiqc.html --title "Read quality" --force
```

Запустим скрипт

```bash
$ bash scripts/qs_script.sh
```

- *Результаты*
    
    [multiqc.html](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/multiqc.html)
    

# Картирование

## Скачаем файлы с референсом

Составим скрипт для скачивания

```bash
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
```

Запустим скрипт

```bash
$ bash scripts/load_reference.sh
```

## Hisat2

Скачаем `hisat2` в наше окружение

```bash
$ conda install hisat2 -c bioconda
```

Проиндексируем геном

```bash
$ mkdir -p references/hisatIndex # создадим папку для индексов
$ hisat2-build -p 8 references/genome.fa references/hisatIndex/musgrcm39Index # запустим индексирование
```

Теперь запустим картирование образцов

```bash
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
```

Запустим скрипт

```bash
$ bash scripts/run_hisat_alignment.sh
```

Для дальнейшей работы нам нужно перевести файлы из формата `.sam` в формат `.bam`, отсортировать их и проиндексировать. Чтобы это сделать воспользуемся `samtools`.

```bash
$ conda install samtools
```

<aside>
🖊️ Не будем выбрасывать риды, которые откартировались неуникально, т.к. `htseq` самостоятельно не будет учитывать их при подсчёте экспрессии (только если его не попросить обратное)

</aside>

- *P.S. как* `htseq` *считает риды*
    
    Документация: [ссылка](https://htseq.readthedocs.io/en/latest/htseqcount.html)
    
    ![Untitled](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled.png)
    

Напишем скрипт. Положим отсортированые `.bam` файлы и их индексы в отдельную папку файлы в отдельную папку.

```bash
#!/bin/bash

mkdir -p sortedbam

for infile in hisatOutput/*
do
  sample=$(basename -a ${infile})
  samtools view -@4 -b hisatOutput/${sample}/alignment.sam > hisatOutput/${sample}/alignment.bam
  samtools sort -@4 hisatOutput/${sample}/alignment.bam -o sortedbam/${sample}_sorted.bam
  samtools index -@4 sortedbam/${sample}_sorted.bam
done
```

Запустим скрипт

```bash
$ bash scripts/sambam_sort_index.sh
```

## Kallisto

В качестве дополнительного задания сделаем квантификацию экспрессии при помощи `kallisto`, который выполняет псевдовыравнивание.

Скачаем [kallisto](http://pachterlab.github.io/kallisto/) в наше окружение. Я здесь скачиваю версию 0.48.0 потому-что более новые версии у меня не работают. Сталкиваюсь со следующей ошибкой https://github.com/pachterlab/kallisto/issues/399

```bash
$ conda install -c conda-forge -c bioconda kallisto==0.48.0
```

Теперь нам нужно проиндексировать транскриптом. 

```bash
$ kallisto index -i references/musGRCm39_kallisto references/transcriptome.fa
```

Теперь напишем скрипт для квантификации наших ридов при помощи kallisto. Я не нашёл информацию о длине фрагментов в описании эксперимента на [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80550). Поэтому будем запускать со стандартными параметрами. Будем сохранять логи в отдельный файл чтобы вклюить их в отчёт `multiqc`

```bash
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
```

<aside>
🖊️ For single end reads, you **must** supply the length and standard deviation of the fragment length (not the read length) ([quick note about the difference between fragment length and read length](https://www.ecseq.com/support/ngs/why-do-the-reads-all-have-the-same-length-when-sequencing-differently-sized-fragments)). The exact information can be obtained by Bioanalyzer/Fragment Analyzer results on the prepared RNA-seq libraries. If you do not have this information, using fragment length = 180 and sd = 20 (like Kallisto manual page) could be a good start.

</aside>

Запустим скрипт

```python
$ bash scripts/kallisto_quantify.sh
```

- Результаты
    
    [multiqcfinal.html](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/multiqcfinal.html)
    

# Alignment QC

Теперь давайте проверим качество выравнивания. Создадим новое окружение для RSeQC

```bash
$ conda create -n rseqc
$ conda activate rseqc
$ conda install -c bioconda rseqc
```

Для запуска самой интересной функции `geneBody_coverage` нам понадобится сделать кастомный `.bed` файл на основе нашего `.gtf` файла. Более того, `geneBody_coverage` работает очень долго (т.к. написана на `python`), поэтому мы бы хотели запустить её на небольшой подвыборке генов. Будем запускать на генах домашнего хозяйства. Для этого возьмём их Refseq идентификаторы из файла ([mm10.HouseKeepingGenes.bed.gz](https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10.HouseKeepingGenes.bed.gz/download)) с сайта документации [RseQC](https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/).

<aside>
🖊️ По-хорошему нужно было бы сделать так же как и они:
1. Взять гены домашнего хозяйства для человека
2. Найти для них гомологи у мыши

</aside>

Напишем скрипт чтобы вытащить нужные нам гены из `.bed` файла, который мы скачали, а так же составить наш `.bed` файл на основе [Biomart](https://www.ensembl.org/biomart/martview/87c698ebbceb6b8e05466761bcb6fa0c) аннотации. 

```python
import pandas as pd

# load a list of house-keeping genes
hk_df = pd.read_csv('references/mm10.HouseKeepingGene.bed', sep='\t', header=None)
hk_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'ItemRgb', 'blockCount', 'blockSizes', 'blockStarts']
hk_df = hk_df.set_index('name')
refseq_id = list(hk_df.index.unique())

# load gtf file to get chr positions
gtf  = pd.read_csv("references/genome.gtf", sep='\t', header=None, skiprows=5, usecols=[0, 3, 4, 8], dtype='str')
gtf.columns = ['chr_ind', 'chrstart', 'chrend', 'Gene stable ID']
gtf['Gene stable ID'] = gtf['Gene stable ID'].map(lambda x: x.split(';')[0].replace('gene_id "', '').replace('"', '') if 'transcript_id' not in x else None)
gtf.dropna(inplace=True)

# load biomart annotation to get refseq ids
biomart_df = pd.read_csv("references/biomart_export.txt", sep='\t')
biomart_hk = biomart_df[biomart_df["RefSeq mRNA ID"].isin(refseq_id)]
biomart_hk = biomart_hk.merge(gtf, on='Gene stable ID').set_index('RefSeq mRNA ID')

# create a bed file
merge_col = ['Gene stable ID', 'chrstart', 'chrend', 'chr_ind']
merge_df = hk_df.merge(biomart_hk[merge_col], left_index=True, right_index=True).reset_index()
bed_file = merge_df[['chr_ind', 'chrstart', 'chrend', 'Gene stable ID', 'score', 'strand', 'thickStart', 'thickEnd', 'ItemRgb', 'blockCount', 'blockSizes', 'blockStarts']]
bed_file.to_csv('references/hk_genes_ensembl.bed', sep='\t', header=None, index=None)

```

Запустим скрипт:

```bash
$ bash scripts/create_bed.py
```

- *Результат*
    
    [hk_genes_ensembl.bed](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/hk_genes_ensembl.bed)
    

Теперь запустим выполнение функции `geneBody_coverage`, используя наш кастомный  `.bed` файл для генов домашнего хозяйства.

```bash
$ mkdir -p RSeQCresults/geneBody_coverage
$ geneBody_coverage.py -r references/hk_genes_ensembl.bed -i sortedbam -o RSeQCresults/geneBody_coverage/
```

- *Результат*
    
    [.geneBodyCoverage.curves.pdf](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/.geneBodyCoverage.curves.pdf)
    

Остальные команды будем запускать c помощью скрипта используя `.bed` файл всего генома.

```bash
#!/bin/bash

# make directories for each tool
mkdir -p RSeQCresults/read_distribution
mkdir -p RSeQCresults/infer_experiment
mkdir -p RSeQCresults/junction_saturation
mkdir -p RSeQCresults/junction_annotation

for infile in $(ls hisatOutput/*/*sorted.bam)
do
  sample=$(basename -a -s _sorted.bam ${infile})
  echo "running read distribution on $sample"
  read_distribution.py  -i $infile -r reference/hg38_RefSeq.bed > RSeQCresults/read_distribution/$sample.readCoverage.txt
  infer_experiment.py  -i $infile -r reference/hg38_RefSeq.bed > RSeQCresults/infer_experiment/$sample.infer_experiment.txt
  junction_saturation.py -i $infile -r reference/hg38_RefSeq.bed -o RSeQCresults/junction_saturation/$sample
  junction_annotation.py -i $infile -r reference/hg38_RefSeq.bed -o RSeQCresults/junction_annotation/$sample
done
```

Запустим скрипт

```python
$ bash rseqc_script.sh
```

- Результат
    
    Некоторые профили не подтягиваются в multiqc, поэтому я нарисовал их при помощи питона
    
    ![clipping profile](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%201.png)
    
    clipping profile
    
    ![delition progile](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%202.png)
    
    delition progile
    
    ![insertion profile](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%203.png)
    
    insertion profile
    
    ![duplication profile](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%204.png)
    
    duplication profile
    
    ![mismatch profile](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%205.png)
    
    mismatch profile
    

Запутим `multiqc` чтобы собрать результаты вместе

```bash
$ multiqc . -f -d -dd 2 --outdir qc --filename multirseqc.html
```

- *Результат*
    
    [multirseqc.html](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/multirseqc.html)
    

# Визуализация

## Genome Browser

Заходим на сайт [Genome Browser](https://genome.ucsc.edu/cgi-bin/hgCustom), выбираем в качестве организма мышку и добавляем кастомный трек. По ключу `bigDataUrl` мы передаём ссылку на `.bam` файл, по ключу `bigDataIndex` ссылку на `.bai` файл. Файлы `.bam` и `.bai` я предварительно загрузил на [figshare](https://figshare.com/articles/dataset/Untitled_Item/25481203?file=45279691).

```
track type=bam name="SRR3414630" bigDataUrl=https://figshare.com/ndownloader/files/45279703 bigDataIndex=https://figshare.com/ndownloader/files/45279691
```

- Выглядит оно как-то так
    
    ![Untitled](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%206.png)
    

Получилось что-то такое. Но открывается не для всех хромосом (видимо потому что нужно было выравнивать на геном от ucsc).

![Untitled](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%207.png)

## Persephone

Я загрузил genome.fa, genome.gtf и свой .bam  файл на сайте [Persephone](https://web.persephonesoft.com/). И получилось что-то такое. Всё делается очень просто и легко, советую.

![Untitled](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%208.png)

![Untitled](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/Untitled%209.png)

# Квантификация

Установим пакет используя конду. Для этого создадим новое окружение для него

```bash
$ conda create -n htseq
$ conda acivate htseq
$ conda install HTSeq
```

Напишем скрипт для подсчёта экспрессии при помощи `htseq-count`

```bash
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
```

Запустим скрипт

```python
$ bash scripts/htseq_quantify.sh
```

Запутим `multiqc` чтобы собрать результаты вместе

```bash
$ multiqc . -f -d -dd 2 --outdir qc --filename multiqcfinal.html
```

- *Результат*
    
    [multiqcfinal.html](RNA-seq%202cabfaf5caa9455f98df9839ec0ceb01/multiqcfinal.html)
    

# Дифференциальная экспрессия

## Нормализация

## Сравнение экспрессий