#!/bin/bash

#####################
# READ BEFORE RUNNING
#
# this shell script runs the covid-seq pipeline using slurm.
# COPY script to the directory you want to run it in.
# make sure raw reads are in a directory named reads
# edit covid reference DATABASE name below, this should be located in /gs1/RTS/NextGen/ngs_dbs/virus/Coronavirus
####################

############################################
#change this for to change the name of the output file

OUTPUT="NGS_05XX_results.xlsx"

############################################

ALIGNER="bowtie2"
#ALIGNER="bwa"

############################################
#for NextSeq runs use min length of 50

LENGTH=75
#LENGTH=50

############################################

PLOIDY=2
#PLOIDY=5

############################################

#SRUN="srun -p bigmem"
#SRUN="srun -p bigmem,int"
SRUN="srun -p int"

############################################

THREADS=8

############################################
# Select "lig" script if Illumina RNA-seq ligation kit was used

#VCPIPE="covid_pipeline.v6.pl"
VCPIPE="covid_pipeline_lig.v6.pl"

############################################
#pick database to use based on Variant Classification Request


## Alpha (B.1.1.7 and Q lineages; earliest sequences from UK)
DATABASE="WA1"             #19B
#DATABASE="APU-POLBA01"
#DATABASE="RML_18497_UK"
#DATABASE="MA10"            #19A

## Beta (B.1.351 and descendent lineages; earliest sequences from South Africa)
#DATABASE="RML_18500_SA"
#DATABASE="RML_18867_MD"

## Delta (B.1.617.2 and AY lineages)
#DATABASE="EPI_ISL_1823618"  #21J
#DATABASE="RML7"
#DATABASE="MZ082533"

## Epsilon (B.1.427 and B.1.1429)

## Eta (B.1.525)

## Gamma (P.1 and descendent lineages; earliest sequences from Brazil)

## Kappa (B.1.617.1)
#DATABASE="EPI_ISL_1675223"

## 1.617.3

## Mu (B.1.621, B.1.621.1)

## Zeta (P.2; earliest sequences from Brazil)

## Omicron (B.1.1.529)
#DATABASE="EPI_ISL_6939036"
## Omicron BA.2
#DATABASE="BA2"
## Omicron BA.5
#DATABASE="SCV2_USA_COR-22-063113_2022_BA.5"


#DATABASE="USA-WA1"      #annotations are orf; suggest to use WA1
#DATABASE="CA_SEARCH"    #contains Ns don't use
#DATABASE="Eng"          #contains Ns don't use
#DATABASE="Eng2"         #contains Ns don't use
#DATABASE="SA_KRISP"     #contains Ns don't use
#DATABASE="MD-HP01541"   #contains Ns don't use

#############################################


mkdir -p consensus
# List your input files

READS=`ls reads/*_R1_001.fastq`

export TMPDIR=~/tmp/


for READ in $READS
do


	SNUM=`echo $READ | grep -P -o '\w+_R'`

	echo $SNUM


	$SRUN -c $THREADS $VCPIPE -f ${SNUM}1_001.fastq -r ${SNUM}2_001.fastq -n $THREADS -a $ALIGNER -d $DATABASE -p $PLOIDY -l $LENGTH &
done
wait

# Wait for all jobs to finish before exiting.

$SRUN bash -c "cat */*_vir_counts.txt > All_vir_hits.txt"
$SRUN bash -c "cat */*_prok_counts.txt > All_prok_hits.txt"
$SRUN bash -c "cat */*_covid_counts.txt > All_covid_hits.txt"

$SRUN bash -c "grep -c '@M0' reads/*R1_001.fastq >> read_counts.txt"
$SRUN bash -c "grep -c '@M0' reads/*sorted_for.fastq >> read_counts.txt"

mkdir -p vcf
cp 1*/*vcf* vcf/
cp 2*/*vcf* vcf/
cp LIB*/*vcf* vcf/



VCFS=`ls vcf/*filtered.vcf`


for VCF in $VCFS
do


	SNUM=`echo $VCF | grep -P -o '\w+.filtered.vcf'`

	echo $SNUM


	VARIANT+=" --variant vcf/$SNUM"

done
wait

echo $VARIANT


$SRUN -c 2 gatk3 -T CombineVariants -R /gs1/RTS/NextGen/ngs_dbs/virus/Coronavirus/${DATABASE}.fa $VARIANT -o all.filtered.vcf

java -jar /gs1/RTS/NextGen/snpEff/snpEff.jar -no-intergenic -no-upstream -no-downstream -no-utr $DATABASE all.filtered.vcf > all.filtered.ann.vcf

$SRUN  -c 2 gatk VariantsToTable -F CHROM -F POS -F TYPE -F REF -F ALT -F ANN -GF AD -V all.filtered.ann.vcf -O all_filtered_variants_ann.txt

$SRUN reformat_snpEff_results.pl -i all_filtered_variants_ann.txt -o final_filtered_variants.txt




RAWS=`ls vcf/*raw.vcf`



for RAW in $RAWS
do


	SNUM=`echo $RAW | grep -P -o '\w+.raw.vcf'`

	echo $SNUM


	VARS+=" --variant vcf/$SNUM"

done
wait

echo $VARS


$SRUN -c 2 gatk3 -T CombineVariants -R /gs1/RTS/NextGen/ngs_dbs/virus/Coronavirus/${DATABASE}.fa $VARS -o all.raw.vcf

java -jar /gs1/RTS/NextGen/snpEff/snpEff.jar -no-intergenic -no-upstream -no-downstream -no-utr $DATABASE all.raw.vcf > all.raw.ann.vcf

$SRUN -c 2 gatk VariantsToTable -F CHROM -F POS -F TYPE -F REF -F ALT -F ANN -GF AD -V all.raw.ann.vcf -O all_raw_variants_ann.txt

$SRUN reformat_snpEff_results.pl -i all_raw_variants_ann.txt -o final_raw_variants.txt


$SRUN covid2excel.pl -o $OUTPUT

#if you want to rename consensus seq in fastq file, you need to make 2 changes to below:
# change SVG_040 to name you want to use and change SampleID
##### rename fasta file
#perl -p -i -e 's/>1/>SVG-040/' consensus/20894_consensus.fa
#perl -p -i -e 's/>1/>SVG-041/' consensus/20895_consensus.fa
#perl -p -i -e 's/>1/>SVG-042/' consensus/20896_consensus.fa
#perl -p -i -e 's/>1/>SVG-043/' consensus/20897_consensus.fa

