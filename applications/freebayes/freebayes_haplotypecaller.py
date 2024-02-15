# ################################### INFO ######################################## #
# FASTP: trim fastq/fasta
# https://github.com/OpenGene/fastp
# ################################### IMPORT ###################################### #

import os
import sys
import json
import glob
import collections
from inspect import currentframe, getframeinfo
# ################################### INCLUDE ##################################### #
# ################################### FUNCTIONS ################################### #


def freebayes_haplotypecaller_script(general_Dict):
	"""
	https://github.com/freebayes/freebayes
	"""
	IO_Dict = {}
	execution_script = ""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	PLATFORM = general_Dict["CONFIG"]["EXECUTION_PLATFORM"]
	#TARGET_REFERENCE_DICT = general_Dict["CONFIG"]["REFERENCE"]["TARGET_REFERENCE"][PLATFORM]
	TARGET_REFERENCE_LIST = general_Dict["CONFIG"]["TARGET_LIST"]
	TARGET_REFERENCE_PATH = general_Dict["CONFIG"]["TARGET_REFERENCE_PATH"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# ++++++++++++++++++++++++++++++++++++++++++++\n"
	if PLATFORM == "BIOWULF":
		execution_script += "module load freebayes\n"
		execution_script += "module load bcftools\n"
		execution_script += "module load rtg-tools\n"
		execution_script += "module load GATK/4.2.5.0\n"
		execution_script += "module load snpEff\n"
		execution_script += "FREEBAYES=$(which freebayes)\n"
		execution_script += "BCFTOOLS=$(which bcftools)\n"
		execution_script += "GATK=$(which gatk)\n"
		execution_script += "BGZIP=$(which bgzip)\n"
		execution_script += "RTG=$(which rtg)\n"
		execution_script += "SNPEFF=\"java -Xmx12g -jar $SNPEFF_JAR\"\n"
	elif PLATFORM == "BIGSKY":
		execution_script += "FREEBAYES=$(which freebayes)\n"
		execution_script += "BCFTOOLS=$(which bcftools)\n"
		execution_script += "BGZIP=$(which bgzip)\n"
		execution_script += "RTG=/gs1/RTS/NextGen/bin/grs_virmap/bin/rtg\n"
		execution_script += "GATK=$(which gatk)\n"
		execution_script += "SNPEFF=\"java -jar /gs1/RTS/NextGen/snpEff/snpEff.jar\"\n"
	execution_script += "# ++++++++++++++++++++++++++++++++++++++++++++\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for each_input in general_Dict["INPUT"]:
		##
		IO_Dict["INPUT1"] = each_input
		PREFIX = os.path.basename(each_input).replace(".bam", "")
		TARGET = PREFIX.split("bowtie2_map.")[1]
		IO_Dict["OUTPUT1"] = general_Dict["OUTPUT_PATH"] + PREFIX + ".freebayes_haplotypecaller.vcf.gz"
		IO_Dict["OUTPUT2"] = general_Dict["OUTPUT_PATH"] + PREFIX + ".freebayes_haplotypecaller.txt"
		IO_Dict["OUTPUT3"] = general_Dict["OUTPUT_PATH"]  + PREFIX + ".freebayes_haplotypecaller.snpeff.vcf.gz"
		IO_Dict["OUTPUT4"] = general_Dict["OUTPUT_PATH"] + PREFIX + ".freebayes_haplotypecaller.snpeff.txt"
		#+
		IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + PREFIX + ".freebayes_haplotypecaller.bcftools_stats.txt"
		IO_Dict["REPORT2"] = general_Dict["REPORT_PATH"] + PREFIX + ".freebayes_haplotypecaller.rtg_vcfstats.txt"
		IO_Dict["REPORT3"] = general_Dict["REPORT_PATH"] + PREFIX + ".freebayes_haplotypecaller.snpeff_summary.csv"
		IO_Dict["REPORT4"] = general_Dict["REPORT_PATH"] + PREFIX + ".freebayes_haplotypecaller.snpeff_summary.html"
		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OUTPUT = general_Dict["CONFIG"]["OUTPUT"]
		TITLE = general_Dict["CONFIG"]["TITLE"]
		IO_Dict["TARGET_REF_NAME"] = TARGET
		IO_Dict["TARGET_REF"] = TARGET_REFERENCE_PATH + "/" + TARGET + "/" + TARGET + ".fa"
		#IO_Dict["TARGET_REF"] = TARGET_REFERENCE_DICT[TARGET].split(".fa")[0] + ".fa"
		IO_Dict["OUTPUT_PATH"] = general_Dict["OUTPUT_PATH"]
		IO_Dict["SNPEFF_CONFIG"] = TARGET_REFERENCE_PATH + "/snpEff.config"
		IO_Dict["NCORE"] = 5
		IO_Dict["freebayes_haplotypecaller_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["VARIANT_CALLING"]["freebayes_haplotypecaller"]
		IO_Dict["snpeff_ann_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["VARIANT_CALLING"]["snpeff_ann"]
		execution_script += "# +++++++++++++++++++++++++++++++++++++++++++++\n"
		execution_script +="""
			$FREEBAYES \\
			-p {NCORE} \\
			{freebayes_haplotypecaller_parameters} \\
			{INPUT1} \\
			-f {TARGET_REF} \\
			> {OUTPUT1}.tmp \\
			2>> {LOG_FILE}

			$BGZIP \\
			-c {OUTPUT1}.tmp \\
			> {OUTPUT1} \\
			2>> {LOG_FILE}
			
			$BCFTOOLS index -f  \\
			{OUTPUT1} >> \\
			{LOG_FILE} 2>&1

			$GATK IndexFeatureFile \\
			-I {OUTPUT1} \\
			>> {LOG_FILE} 2>&1

			$GATK VariantsToTable \\
			-F CHROM -F POS -F TYPE -F REF -F ALT -GF AD \\
			-R {TARGET_REF} \\
			-V {OUTPUT1} \\
			-O {OUTPUT2} \\
			>> {LOG_FILE} 2>&1

			# $BCFTOOLS query --print-header \\
			# -f \"%CHROM\\t%POS\\t%TYPE\\t%REF\\t%ALT[\\t%GT\\t%RO\\t%AO]\" \\
			# {OUTPUT1} > {OUTPUT2} 2>> {LOG_FILE}


			$BCFTOOLS stats {OUTPUT1} \\
			> {REPORT1} 2>> {LOG_FILE}

			$RTG vcfstats {OUTPUT1} \\
			> {REPORT2} 2>> {LOG_FILE}

			$SNPEFF ann \\
			{snpeff_ann_parameters} \\
			-config {SNPEFF_CONFIG} -v {TARGET_REF_NAME} -csvStats {REPORT3} -stats {REPORT4}  \\
			{OUTPUT1} > {OUTPUT3}.tmp 2>> {LOG_FILE} 

			$GATK VariantsToTable \\
			-F CHROM -F POS -F TYPE -F REF -F ALT -F ANN -GF AD  \\
			-R {TARGET_REF} \\
			-V {OUTPUT3}.tmp \\
			-O {OUTPUT4} \\
			>> {LOG_FILE} 2>&1

			$BGZIP \\
			-c {OUTPUT3}.tmp \\
			> {OUTPUT3} \\
			2>> {LOG_FILE}


			$BCFTOOLS index -f \\
			{OUTPUT3} >> \\
			{LOG_FILE} 2>&1
			

			
			#rm {OUTPUT_PATH}*.tmp
		""".format(**IO_Dict)
		execution_script += "# --------------------------------------------\n"

	else:
		##
		pass
		
	return execution_script


def build_freebayes_haplotypecaller_execution_script(general_Dict):
	"""
	"""
	execution_script = freebayes_haplotypecaller_script(general_Dict)
	return execution_script

# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #