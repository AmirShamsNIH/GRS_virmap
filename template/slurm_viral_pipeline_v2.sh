#!/bin/bash

#####################
# READ BEFORE RUNNING
#
# this shell script runs the vircapseq pipeline using slurm.
# COPY script to the directory you want to run it in.
# make sure raw reads are in a directory named reads
# edit vircapseq database name below, this should be located in /gs1/RTS/NextGen/ngs_dbs/virus/VirCapSeq
####################

############################################
#change this for to change the name of the output file

OUTPUT="GRS_0167_deWit_results.xlsx"
#OUTPUT="xxxxxxxxxxxxxxx.xlsx"

############################################

REFERENCE="/gs1/RTS/NextGen/ngs_dbs/virus/Nipah/Nipah_virus_Bangladesh_LVMiSeqConsensus_20181119.fa"
DATABASE="Nipah_virus_Bangladesh_LVMiSeqConsensus_20181119"

############################################

ALIGNER="bowtie2"
#ALIGNER="bwa"

############################################
#for NextSeq runs use min length of 50

LENGTH=75
#LENGTH=50

############################################

#SRUN="srun -p bigmem"
#SRUN="srun -p bigmem,int"
SRUN="srun -p int"

############################################

THREADS=8
#VCPIPE="viral_pipeline.pl "
VCPIPE="viral_pipeline_lig.pl "

############################################


# List your input files

READS=`ls reads/*_R1_001.fastq`

export TMPDIR=~/tmp/




for READ in $READS
do


	SNUM=`echo $READ | grep -P -o '\w+_R'`
	
	echo $SNUM


	$SRUN -c $THREADS $VCPIPE -f ${SNUM}1_001.fastq -r ${SNUM}2_001.fastq -d $REFERENCE -n $THREADS -a $ALIGNER -l $LENGTH &
done
wait

# Wait for all jobs to finish before exiting.

$SRUN bash -c "cat */*_vir_counts.txt > All_vir_hits.txt"
$SRUN bash -c "cat */*_prok_counts.txt > All_prok_hits.txt"
$SRUN bash -c "cat */*_viral_counts.txt > All_viral_hits.txt"
$SRUN bash -c "cat */*_viral_counts.txt > All_covid_hits.txt"

$SRUN bash -c "grep -c '@M0' reads/*R1_001.fastq >> read_counts.txt"
$SRUN bash -c "grep -c '@M0' reads/*sorted_for.fastq >> read_counts.txt"

mkdir vcf
#cp 1*/*vcf* vcf/
#cp 2*/*vcf* vcf/
cp LIB*/*vcf* vcf/



VCFS=`ls vcf/*filtered.vcf.gz`


for VCF in $VCFS
do


	SNUM=`echo $VCF | grep -P -o '\w+.filtered.vcf.gz'`
	
	echo $SNUM


	VARIANT+=" vcf/$SNUM"

done
wait

echo $VARIANT

$SRUN  -c 4 bcftools merge --threads 4 -m none -O v -o all.filtered.vcf $VARIANT

java -jar /gs1/RTS/NextGen/snpEff/snpEff.jar -no-intergenic -no-upstream -no-downstream -no-utr $DATABASE all.filtered.vcf > all.filtered.ann.vcf

$SRUN  -c 2 gatk VariantsToTable -F CHROM -F POS -F TYPE -F REF -F ALT -F ANN -GF AD -V all.filtered.ann.vcf -O all_filtered_variants_ann.txt

$SRUN reformat_snpEff_results.pl -i all_filtered_variants_ann.txt -o final_filtered_variants.txt

$SRUN -c 4 gatk VariantsToTable -F CHROM -F POS -F TYPE -F REF -F ALT -GF AD -V all.filtered.vcf -O all_filtered_variants.txt


RAWS=`ls vcf/*raw.vcf.gz`

for RAW in $RAWS
do


	SNUM=`echo $RAW | grep -P -o '\w+.raw.vcf.gz'`
	
	echo $SNUM


	VARS+=" vcf/$SNUM"

done
wait

echo $VARS


$SRUN  -c 4 bcftools merge --threads 4 -m none -O v -o all.raw.vcf $VARS

java -jar /gs1/RTS/NextGen/snpEff/snpEff.jar -no-intergenic -no-upstream -no-downstream -no-utr $DATABASE all.raw.vcf > all.raw.ann.vcf

$SRUN -c 2 gatk VariantsToTable -F CHROM -F POS -F TYPE -F REF -F ALT -F ANN -GF AD -V all.raw.ann.vcf -O all_raw_variants_ann.txt

$SRUN reformat_snpEff_results.pl -i all_raw_variants_ann.txt -o final_raw_variants.txt

$SRUN  -c 4 gatk VariantsToTable -F CHROM -F POS -F TYPE -F REF -F ALT -GF AD -V all.raw.vcf -O all_raw_variants.txt

#$SRUN viral2excel.pl -o $OUTPUT

$SRUN covid2excel.pl -o $OUTPUT

