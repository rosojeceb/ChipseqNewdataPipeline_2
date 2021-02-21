### Authors: Francisco Manuel Gordillo Cantón, Marta Pérez de Miguel, Rosario Ojeda Ceballos
###
### Contact: fcomanuelgordilloc@gmail.com, mrtprz@gmail.com, rosojeceb@gmail.com
###
### ChIP-Seq Data Analysis
###
### Pipeline Script

#! /bin/bash

## Firstly, a help message appears when no parameter is provided

if [ $# -ne 1 ]
then
 echo "Usage: bash ChipseqNewdataPipeline.sh  <params_file>"
 echo ""
 echo "params.file: Input file with the specific parameters. Please, check the file test_params.txt to visualize an example"
 exit 
fi

##Aquí creo que podríamos o describir los parámetros o sugerir que el usuario vaya a params para verlo ahí

## Secondly, the relevant parameters are introduced

PARAMS=$1

echo "======================"
echo "| LOADING PARAMETERS |"
echo "======================"

INSDIR=$(grep installation_directory: $PARAMS | awk '{ print $2 }')
echo "Installation directory = $INSDIR"

WD=$(grep working_directory: $PARAMS | awk '{ print $2 }')
echo "Working directory = $WD"

EXP=$(grep experiment_name: $PARAMS | awk '{ print $2 }')
echo "Experiment name = $EXP"

NUMSAMPLES=$(grep number_samples: $PARAMS | awk '{ print $2 }')
echo "Number of samples = $NUMSAMPLES"

CHIP_NUMSAMPLES=$(grep number_chip_samples: $PARAMS | awk '{ print $2 }' )
echo "Number of chip samples = $CHIP_NUMSAMPLES"

INPUT_NUMSAMPLES=$(grep number_input_samples: $PARAMS | awk '{ print $2 }' )
echo "Number of input samples = $INPUT_NUMSAMPLES"

GENOME=$(grep path_genome: $PARAMS | awk '{ print $2 }')
echo "Reference genome = $GENOME"

ANNOT=$(grep path_annotation: $PARAMS | awk '{ print $2 }')
echo "Annotation of the genome = $ANNOT"

PROM=$(grep promoter_length: $PARAMS | awk '{ print $2 }' )


CHIP_SAMPLES=()
i=0
while [ $i -lt  $CHIP_NUMSAMPLES ]
do
  j=$(( i + 1 ))
  CHIP_SAMPLES[$i]=$(grep path_chip_sample_$j: $PARAMS | awk '{ print $2 }')
  ((i++))
done
echo "Chip samples = "
echo ${CHIP_SAMPLES[@]}

INPUT_SAMPLES=()
i=0
while [ $i -lt  $INPUT_NUMSAMPLES ]
do
  j=$(( i + 1 ))
  INPUT_SAMPLES[$i]=$(grep path_input_sample_$j: $PARAMS | awk '{ print $2 }')
  ((i++))
done
echo "Input samples = "
echo ${INPUT_SAMPLES[@]}

## Thirdly, a work space is generated

echo""
echo "=========================="
echo "| CREATING WORKING SPACE |"
echo "=========================="
echo""

cd $WD
mkdir $EXP
cd $EXP
mkdir genome annotation results samples scripts

# Obtaining the reference genome  

cd genome
cp $GENOME genome.fa

# Obtaining the annotation

cd ../annotation
cp $ANNOT annotation.gtf


# Creating reference genome index
cd ../genome
bowtie2-build genome.fa index
echo "Genome index created!"

cd ../samples
mkdir chip input

cd chip

i=1
while [ $i -le  $CHIP_NUMSAMPLES ]
do
  j=$(($i - 1))
  cp ${CHIP_SAMPLES[$j]} chip_sample_$i.fq.gz
  ((i++))
done

cd ../input

i=1
while [ $i -le  $INPUT_NUMSAMPLES ]
do
  j=$(($i - 1))
  cp ${INPUT_SAMPLES[$j]} input_sample_$i.fq.gz
  ((i++))
done


cd ../../results
echo "Processing individual samples..."
i=1
while [ $i -le $CHIP_NUMSAMPLES ]
do
  bash $INSDIR/chip_sample_processing.sh $WD/$EXP/samples/chip $i $CHIP_NUMSAMPLES $INSDIR $NUMSAMPLES $WD $EXP $PROM
  ((i++))
done

i=1
while [ $i -le $INPUT_NUMSAMPLES ]
do
  bash $INSDIR/input_sample_processing.sh $WD/$EXP/samples/input $i $INPUT_NUMSAMPLES $INSDIR $NUMSAMPLES $CHIP_NUMSAMPLES $WD $EXP $PROM
  ((i++))
done


cd $WD/$EXP/results

#Finding Peaks with macs2 function
echo "Finding peaks..."

i=1
while [ $i -le $CHIP_NUMSAMPLES ]
do
   mkdir peaks_$i
   cd peaks_$i
   macs2 callpeak -t ../samples/chip/chip_sample_$i.bam -c ../samples/input/input_sample_$i.bam -f BAM --outdir . -n PEAK_$i
   cd ..
   ((i++))
done

echo ""
echo "All peaks found!"
echo ""

##Con --known model se podría hacer para marcas epigenéticas (no genera modelo en la llamada a picos)
##macs2 genera todos los archivos con nombre PEAK_i (no sé qué nombre poner a los picos si no)

#HOMER analysis
cd $WD/$EXP/results
mkdir homer
cd homer
i=1
while [ $i -le $CHIP_NUMSAMPLES ]
do
   mkdir motifs_$i
   cd motifs_$i
   #findMotifsGenome.pl $WD/$EXP/results/PEAK_${i}_peaks.narrowPeak tair10 $WD/$EXP/results/homer -size 100 -len 8
   cd ..
   ((i++))
done

echo ""
echo "Peaks analyzed with HOMER"
echo ""

#R data analysis
cd $WD/$EXP/results
mkdir Rresults
cd Rresults
echo "The promoter length is $PROM"
i=1
while [ $i -le $CHIP_NUMSAMPLES ]
do
   mkdir Rresults_$i
   cd Rresults_$i
   Rscript $INSDIR/chipseq_R_analysis.R $PROM $WD/$EXP/results/PEAK_${i}_peaks.narrowPeak $WD/$EXP/results/PEAK_${i}_summits.bed
   cd ..
   ((i++))
done

echo "ChIP analysis is done! Thank you for using our pipeline!"