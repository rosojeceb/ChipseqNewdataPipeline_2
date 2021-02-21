### Authors: Francisco Manuel Gordillo Cantón, Marta Pérez de Miguel, Rosario Ojeda Ceballos
###
### Contact: fcomanuelgordilloc@gmail.com, mrtprz@gmail.com, rosojeceb@gmail.com
###
### ChIP-Seq Data Analysis
###
### Input samples processing script



#! /bin/bash

# Reading the parameters

SAMPLE_DIR=$1
i=$2
INPUT_NUMSAMPLES=$3
INSDIR=$4
NUM_SAMPLES=$5
CHIP_NUMSAMPLES=$6
WD=$7
EXP=$8
PROM=$9


echo ""
echo "Processing sample $i..."
echo ""

cd ${SAMPLE_DIR}

## Sample quality control and read mapping to reference genome
fastqc input_sample_$i.fq.gz
bowtie2 -x ../../genome/index -U input_sample_$i.fq.gz -S input_sample_$i.sam

## Generating sorted bam file
samtools sort -o input_sample_$i.bam input_sample_$i.sam
rm input_sample_$i.sam
samtools index input_sample_$i.bam
