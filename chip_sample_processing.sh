### Authors: Francisco Manuel Gordillo Cantón, Marta Pérez de Miguel, Rosario Ojeda Ceballos
###
### Contact: fcomanuelgordilloc@gmail.com, mrtprz@gmail.com, rosojeceb@gmail.com
###
### ChIP-Seq Data Analysis
###
### Chip samples processing script


#! /bin/bash

# Reading the parameters

SAMPLE_DIR=$1
i=$2
CHIP_NUMSAMPLES=$3
INSDIR=$4
NUM_SAMPLES=$5
WD=$6
EXP=$7
PROM=$8

echo ""
echo "Processing sample $i..."
echo ""

cd ${SAMPLE_DIR}

## Sample quality control and read mapping to reference genome
fastqc chip_sample_$i.fq.gz
bowtie2 -x ../../genome/index -U chip_sample_$i.fq.gz -S chip_sample_$i.sam

## Generating sorted bam file
samtools sort -o chip_sample_$i.bam chip_sample_$i.sam
rm chip_sample_$i.sam
samtools index chip_sample_$i.bam