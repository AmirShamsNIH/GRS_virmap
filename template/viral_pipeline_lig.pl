#!/usr/bin/perl

#############################
#
#  martensc 12/10/2020
#
# this script runs a pipeline to analyze covid/capture data
#
#############################

use Getopt::Long;

#get program name
$numSplit = (split /\//,$0);
$name = (split /\//,$0)[$numSplit - 1];

$usage = "Usage:  $name -f <FORWARD READS> -r <REVERSE READS> -d <REFERENCE> -n <NUM THREADS> -a <ALIGNER> -l <MIN LENGTH> \n        $name -h for help\n";

# Set up the command line options.
my $ret = GetOptions ("r|reverse:s", "f|forward:s", "d|database:s", "l|length:i", "n|threads:i", "a|aligner:s", "h|help!");
my $fname = $opt_f;
my $rname = $opt_r;
my $reffa = $opt_d;
my $minlen = $opt_l;
my $numthreads = $opt_n;
my $aligner = $opt_a;
my $help = $opt_h;

# Print help screen if -h option is selected
if ($help or !$fname or !$rname or !$numthreads or !$aligner or !$reffa or !$minlen) {
        print "\n\tThis script runs the vircapseq pipeline ";
        print "\n";
        print "\n$usage";
        print "\n$name arguments:\n\n";
        print "\t-f - Name of forward reads in fastq.\n";
        print "\t-r - Name of reverse reads in fastq.\n";
        print "\t-d - Name of sequence\n";
        print "\t-n - Number of threads to use";
        print "\t-a - Alignment program to use: bowtie2 or bwa.\n";
        print "\t-h - display this help file\n\n";
        exit;
}


########################
#set up constants

my $read_dir = "reads";

#my $adapt = "AGATCGGAAGAGC";
my $adapt = "CTGTCTCTTATACACATCT";  # illumina ligation adapter

#$minlen = 75;

#$numthreads = 16;

my $humRibodb = "/gs1/RTS/NextGen/ngs_dbs/vertebrate/Homo_sapiens/hg38/bowtie2/humRibosomal";
my $hg38db = "/gs1/RTS/NextGen/ngs_dbs/vertebrate/Homo_sapiens/hg38/hisat2/GRCh38";
my $phix = "/gs1/RTS/NextGen/ngs_dbs/virus/PhiX/phiX";
my $repeat = "/gs1/RTS/NextGen/ngs_dbs/other/repeat";

my $blastdb = "/gs1/RTS/NextGen/blast_dbs/vir/ref_viruses_rml";

$coviddb = "/gs1/RTS/NextGen/ngs_dbs/virus/CCHFV/ma-cchfv";
$covidfa = "/gs1/RTS/NextGen/ngs_dbs/virus/CCHFV/ma-cchfv.fa";


$refdb = $reffa;
$refdb =~ s/\.fasta//;
$refdb =~ s/\.fa//;



$virusdb = "/gs1/RTS/NextGen/ngs_dbs/virus/VirCapSeq/ref_genomes";
$prokdb = "/gs1/RTS/NextGen/ngs_dbs/bacteria/VirCapSeq/ref_prok_genomes";

########################

#check to see if files are same sampleID

@forname = split(/\_S/,$fname);
@revname = split(/\_S/,$rname);

if ($forname[0] != $revname[0]) {
	print "reads don't match\n";
	exit;
}

########################

#set file names

$basename = $forname[0];
$output_dir = $basename;
`mkdir -p $output_dir`;


#fastq

$for_in = $read_dir."/".$fname;
$rev_in = $read_dir."/".$rname;
$for_trimmed = $read_dir."/".$basename."_R1_trim.fastq";
$rev_trimmed = $read_dir."/".$basename."_R2_trim.fastq";
$for_trimmed1 = $read_dir."/".$basename."_R1_trim1.fastq";
$rev_trimmed1 = $read_dir."/".$basename."_R2_trim1.fastq";
$for_trimmed2 = $read_dir."/".$basename."_R1_trim2.fastq";
$rev_trimmed2 = $read_dir."/".$basename."_R2_trim2.fastq";
$for_filtered = $read_dir."/".$basename."_R1_trim2_filter.fastq";
$rev_filtered = $read_dir."/".$basename."_R2_trim2_filter.fastq";
$for_filtered2 = $read_dir."/".$basename."_R1_trim2_filter2.fastq";
$rev_filtered2 = $read_dir."/".$basename."_R2_trim2_filter2.fastq";
$sorted = $read_dir."/".$basename."_sorted";
$for_sort = $read_dir."/".$basename."_sorted_for.fastq";
$rev_sort = $read_dir."/".$basename."_sorted_rev.fastq";

$noPhiX =  $read_dir."/".$basename."_noPhiX.fastq";
$noPhiX_for =  $read_dir."/".$basename."_noPhiX.1.fastq";
$noPhiX_rev =  $read_dir."/".$basename."_noPhiX.2.fastq";

$noRep =  $read_dir."/".$basename."_noRep.fastq";
$noRep_for =  $read_dir."/".$basename."_noRep.1.fastq";
$noRep_rev =  $read_dir."/".$basename."_noRep.2.fastq";

$noRibo =  $read_dir."/".$basename."_noRibo.fastq";
$noRibo_for =  $read_dir."/".$basename."_noRibo.1.fastq";
$noRibo_rev =  $read_dir."/".$basename."_noRibo.2.fastq";

$noHg38 =  $read_dir."/".$basename."_noHg38.fastq";
$noHg38_for =  $read_dir."/".$basename."_noHg38.1.fastq";
$noHg38_rev =  $read_dir."/".$basename."_noHg38.2.fastq";

$noCovid =  $read_dir."/".$basename."_noCovid.fastq";
$noCovid_for =  $read_dir."/".$basename."_noCovid.1.fastq";
$noCovid_rev =  $read_dir."/".$basename."_noCovid.2.fastq";

#sam/bam/vcf

$ribosam = $output_dir."/".$basename."_ribo.sam";
$hg38sam = $output_dir."/".$basename."_hg38.sam";
$virsam = $output_dir."/".$basename."_vir.sam";
$repsam = $output_dir."/".$basename."_rep.sam";
$phixsam = $output_dir."/".$basename."_phix.sam";
$proksam = $output_dir."/".$basename."_prok.sam";
$covidsam = $output_dir."/".$basename."_covid.sam";

$nocovidsam = $output_dir."/".$basename."_nocovid.sam";

$covidbam = $output_dir."/".$basename.".bam";
$covidsortbam = $output_dir."/".$basename."_sorted.bam";
$covidsortnodupbam = $output_dir."/".$basename."_sorted_nodup.bam";
$covidsortnoduprgbam = $output_dir."/".$basename."_sorted_nodup_rg.bam";

$virbam = $output_dir."/".$basename."_vir.bam";
$virsortbam = $output_dir."/".$basename."_vir_sorted.bam";
$virsortnodupbam = $output_dir."/".$basename."_vir_sorted_nodup.bam";

$vcf = $output_dir."/".$basename.".raw.vcf";
$vcf_filtered = $output_dir."/".$basename.".filtered.vcf";

$vcf_gz = $output_dir."/".$basename.".raw.vcf.gz";
$vcf_filtered_gz = $output_dir."/".$basename.".filtered.vcf.gz";

#txt files
$vir_count = $output_dir."/".$basename."_vir_counts.txt";
$prok_count = $output_dir."/".$basename."_prok_counts.txt";
$covid_count = $output_dir."/".$basename."_viral_counts.txt";
$picard_met = $output_dir."/".$basename."_picard_met.txt";

#########################

#trim/filter reads

`cutadapt -a $adapt -A $adapt --trim-n -m $minlen -o $for_trimmed -p $rev_trimmed $for_in $rev_in`;

`fastx_trimmer -f 2 -i $for_trimmed -o $for_trimmed1`;
`fastx_trimmer -f 2 -i $rev_trimmed -o $rev_trimmed1`;

`fastq_quality_trimmer -i $for_trimmed1 -o $for_trimmed2 -t 20 -l $minlen`;
`fastq_quality_trimmer -i $rev_trimmed1 -o $rev_trimmed2 -t 20 -l $minlen`;

`fastq_quality_filter -i $for_trimmed2 -o $for_filtered -q 20 -p 80`;
`fastq_quality_filter -i $rev_trimmed2 -o $rev_filtered -q 20 -p 80`;

`remove_n_fastq.pl -i $for_filtered -o $for_filtered2`;
`remove_n_fastq.pl -i $rev_filtered -o $rev_filtered2`;

`sort_fastq.pl -f $for_filtered2 -r $rev_filtered2 -o $sorted`;

#print "sort_fastq.pl -f $for_filtered2 -r $rev_filtered2 -o $sorted\n";

#exit;
#######################

#screen reads

#repeats/phix
`bowtie2 --fast -p $numthreads --no-mixed --no-unal -X 1200 --un-conc $noPhiX -x $phix -1 $for_sort -2 $rev_sort -S $phixsam`;
`bowtie2 --fast -p $numthreads --no-mixed --no-unal -X 1200 --un-conc $noRep -x $repeat -1 $noPhiX_for -2 $noPhiX_rev -S $repsam`; 

#human
`bowtie2 --fast -p $numthreads --no-mixed --no-unal -X 1200 --un-conc $noRibo -x $humRibodb -1 $noPhiX_for -2 $noPhiX_rev -S $ribosam`;
`hisat2 -p $numthreads --no-mixed --no-unal --un-conc $noHg38 -x $hg38db -1 $noRibo_for -2 $noRibo_rev -S $hg38sam`;

#######################

#screen against viruses of interest

#`bowtie2 -p $numthreads --no-mixed --no-unal -X 1500 --un-conc $noCovid -x $coviddb  -1 $noHg38_for -2 $noHg38_rev -S $covidsam`;

if ($aligner eq "bwa") {
	`bwa mem -t $numthreads $reffa $noHg38_for $noHg38_rev > $covidsam`;
	`picard ViewSam I=$covidsam ALIGNMENT_STATUS=Unaligned PF_STATUS=All > $nocovidsam`;
	`picard SamToFastq I=$nocovidsam FASTQ=$noCovid_for SECOND_END_FASTQ=$noCovid_rev`;
}
elsif ($aligner eq "bowtie2") {
	`bowtie2 -p $numthreads --local --very-sensitive --no-unal -X 1500 --un-conc $noCovid -x $refdb  -1 $noHg38_for -2 $noHg38_rev -S $covidsam`;
}
else {
	die "can't fine $aligner\n";
}

`count_genome_hits.pl -i $covidsam -o $covid_count`;

`samtools view -b -o $covidbam $covidsam`;
`samtools sort -o $covidsortbam $covidbam`;
`samtools index $covidsortbam`;
`picard MarkDuplicates REMOVE_DUPLICATES=true I=$covidsortbam O=$covidsortnodupbam M=$picard_met`;
`samtools index $covidsortnodupbam`;

`picard AddOrReplaceReadGroups RGID=$basename RGLB=$basename RGSM=$basename RGPL=illumina RGPU=$basename INPUT=$covidsortnodupbam OUTPUT=$covidsortnoduprgbam`;
`samtools index $covidsortnoduprgbam`;
`gatk HaplotypeCaller -ploidy 2 -R $reffa  -I $covidsortnoduprgbam -O $vcf`;

`bcftools filter -i '(TYPE="snp" && QUAL>500 && INFO/DP>20) || (TYPE="indel" && QUAL>1000 && INFO/DP>20)' -O v -o $vcf_filtered $vcf`;

`bgzip -c $vcf > $vcf_gz`;
`bgzip -c $vcf_filtered > $vcf_filtered_gz`;  

`bcftools index $vcf_gz`;
`bcftools index $vcf_filtered_gz`; 

######################


#####################


#screen against viruses of interest

`bowtie2 -p $numthreads --no-mixed --no-unal -X 1500 -x $virusdb  -1 $noCovid_for -2 $noCovid_rev -S $virsam`;

`count_genome_hits.pl -i $virsam -o $vir_count`;

`samtools view -b -o $virbam $virsam`;
`samtools sort -o $virsortbam $virbam`;
`samtools rmdup $virsortbam $virsortnodupbam`;
`samtools index $virsortnodupbam`;

#########################

#screen against bacteria

`bowtie2 -p $numthreads --no-mixed --no-unal -X 1500 -x $prokdb  -1 $noCovid_for -2 $noCovid_rev -S $proksam`;

`count_genome_hits.pl -i $proksam -o $prok_count`;


#########################



#########################

#file clean up

`rm $for_trimmed`;
`rm $rev_trimmed`;
`rm $for_trimmed2`;
`rm $rev_trimmed2`;
`rm $for_filtered`;
`rm $rev_filtered`;
`rm $for_filtered2`;
`rm $rev_filtered2`;

#`rm $ribosam`;
#`rm $hg38sam`;
`rm $repsam`;
`rm $phixsam`;

`rm $noPhiX_for`;
`rm $noPhiX_rev`;
`rm $noRep_for`;
`rm $noRep_rev`;
`rm $noRibo_for`;
`rm $noRibo_rev`;

######################### 
