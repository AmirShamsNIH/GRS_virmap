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

def custom_virmapDB_script(general_Dict):
	"""
	1)bowtie2 index
	2)index fasta
	3)fasta dictionary
	4)snpEff gff database
	"""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	IO_Dict = {}
	EXECUTION_PLATFORM = general_Dict["CONFIG"]["EXECUTION_PLATFORM"]

	IO_Dict["OUTPUT_PATH"] = general_Dict["OUTPUT_PATH"]
	#+
	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script = "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIOWULF":
		execution_script += "module load samtools\n"
		execution_script += "module load bowtie/2-2.5.1\n"
		execution_script += "module load snpEff\n"
		execution_script += "BOWTIE2=$(which bowtie2)\n"
		execution_script += "SAMTOOLS=\"$(which samtools)\"\n"
		execution_script += "SNPEFF=\"java -Xmx12g -jar $SNPEFF_JAR\"\n"
	elif general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIGSKY":
		execution_script += "SAMTOOLS=$(which samtools)\n"
		execution_script += "BOWTIE2=$(which bowtie2)\n"
		execution_script += "SNPEFF=\"java -jar /gs1/RTS/NextGen/snpEff/snpEff.jar\"\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	execution_script += "# -------------------------------------------\n"
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"

	IO_Dict["bowtie2_build_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["BUILD_ENVIRONEMNT"]["bowtie2_build"]
	IO_Dict["snpeff_build_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["BUILD_ENVIRONEMNT"]["snpeff_build"]
	for each_target_index in general_Dict["CONFIG"]["TARGET_REFERENCE_INDEX"]:
		##
		target_fasta_file = general_Dict["CONFIG"]["TARGET_REFERENCE"][EXECUTION_PLATFORM][each_target_index]["fasta"]
		target_gtf_file = general_Dict["CONFIG"]["TARGET_REFERENCE"][EXECUTION_PLATFORM][each_target_index]["gtf"]
		IO_Dict["TARGET_FASTA"] = target_fasta_file
		IO_Dict["TARGET_GTF"] = target_gtf_file
		gtf_file_extension = os.path.basename(target_gtf_file).split(".")[-1]
		if gtf_file_extension in ["gtf"]:
			IO_Dict["GTF_FORMAT"] = "-gtf22"
			IO_Dict["GTF_EXTENSION"] = ".gtf"
		elif gtf_file_extension in ["gff", "gff3"]:
			IO_Dict["GTF_FORMAT"] = "-gff3"
			IO_Dict["GTF_EXTENSION"] = ".gff"
		elif gtf_file_extension in ["gbk"]:
			IO_Dict["GTF_FORMAT"] = "-genbank"
			IO_Dict["GTF_EXTENSION"] = "." + gtf_file_extension
		IO_Dict["TARGET"] = each_target_index
		execution_script += """
			mkdir -p {OUTPUT_PATH}{TARGET}
			cp {TARGET_FASTA} {OUTPUT_PATH}{TARGET}/{TARGET}.fa
			cd {OUTPUT_PATH}{TARGET}/
			$SAMTOOLS faidx {TARGET}.fa \\
			>> {LOG_FILE} 2>&1
			$SAMTOOLS dict {OUTPUT_PATH}{TARGET}/{TARGET}.fa --output {OUTPUT_PATH}{TARGET}/{TARGET}.dict \\
			>> {LOG_FILE} 2>&1
			
			printf "%s\\n" "# ++++++++++++++++++++++++++++" >> {LOG_FILE}
			cp {TARGET_FASTA} {OUTPUT_PATH}{TARGET}/sequences.fa
			cp {TARGET_GTF} {OUTPUT_PATH}{TARGET}/genes{GTF_EXTENSION}
			
			printf "%s\\n" "{TARGET}.genome: {TARGET}" >> {OUTPUT_PATH}/snpEff.config

			$SNPEFF build {snpeff_build_parameters} {GTF_FORMAT} -dataDir {OUTPUT_PATH} -config {OUTPUT_PATH}/snpEff.config -v {TARGET} \\
			> /dev/null 2>> {LOG_FILE}

			printf "%s\\n" "# ++++++++++++++++++++++++++++" >> {LOG_FILE}
			$BOWTIE2-build --threads {NCORE} {bowtie2_build_parameters} -f {OUTPUT_PATH}{TARGET}/{TARGET}.fa {OUTPUT_PATH}{TARGET}/{TARGET} \\
			>> {LOG_FILE} 2>&1
			
		""".format(**IO_Dict)

	else:
		##
		
		pass

	
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# -------------------------------------------\n"
	return execution_script


def build_custom_virmapDB_execution_script(general_Dict):
	"""
	"""
	execution_script = custom_virmapDB_script(general_Dict)

	return execution_script
# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #