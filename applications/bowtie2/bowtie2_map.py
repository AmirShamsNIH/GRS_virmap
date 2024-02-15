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
# ################################### UTILITY FUNCTIONS ########################### #

# ################################### BOWTIE2 MAP ############################### #

def bowtie2_map_script(general_Dict):
	"""
	https://gensoft.pasteur.fr/docs/bowtie2/2.3.5.1/
	"""
	IO_Dict = {}
	execution_script = ""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	TARGET_REFERENCE_LIST = general_Dict["CONFIG"]["TARGET_LIST"]
	TARGET_REFERENCE_PATH = general_Dict["CONFIG"]["TARGET_REFERENCE_PATH"]
	PLATFORM = general_Dict["CONFIG"]["EXECUTION_PLATFORM"]
	SNAKERULE = general_Dict["SNAKERULE"]
	#TARGET_REFERENCE_DICT = general_Dict["CONFIG"]["REFERENCE"]["TARGET_REFERENCE"][PLATFORM]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if PLATFORM == "BIOWULF":
		execution_script += "module load samtools\n"
		execution_script += "module load bowtie/2-2.5.1\n"
		execution_script += "module load picard\n"
		execution_script += "SAMTOOLS=$(which samtools)\n"
		execution_script += "BOWTIE2=$(which bowtie2)\n"
		execution_script += "PICARD=\"java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar\"\n"
	elif PLATFORM == "BIGSKY":
		execution_script += "SAMTOOLS=$(which samtools)\n"
		execution_script += "BOWTIE2=$(which bowtie2)\n"
		execution_script += "PICARD=$(which picard)\n"
	execution_script += "# -------------------------------------------\n"

	extra_script = "#+++++++++++++++++++++++++++++++++++++++++++++\n"
	extra_script += "rm -rf " + general_Dict["OUTPUT_PATH"] + "*.nodup* \n"
	extra_script += "rm -rf " + general_Dict["OUTPUT_PATH"] + "*.tmp* \n"
	extra_script += "# --------------------------------------------\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for each_target in TARGET_REFERENCE_LIST:
		##
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		OUTPUT = general_Dict["CONFIG"]["OUTPUT"]
		TITLE = general_Dict["CONFIG"]["TITLE"]
		IO_Dict["TARGET_REF"] = TARGET_REFERENCE_PATH + "/" + each_target + "/" + each_target
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if len(general_Dict["INPUT"]) > 1:
			IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
			IO_Dict["INPUT2"] = general_Dict["INPUT"][1]
		else:
			IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
		#+
		IO_Dict["OUTPUT1"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." +  each_target + ".bam"
		IO_Dict["OUTPUT2"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_target + "_mapped.R%.fastq.gz"
		IO_Dict["OUTPUT3"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_target + "_unmapped.R%.fastq.gz"
		#+
		IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_target + ".txt"
		IO_Dict["REPORT2"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_target + ".picard_markduplicate.txt"
		#+
		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
		IO_Dict["SAMPLE"] = general_Dict["SAMPLE"]
		IO_Dict["bowtie2_map_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["ALIGNMENT"]["bowtie2_map"] 
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
		if len(general_Dict["INPUT"]) > 1:
			execution_script += """	
			$BOWTIE2 \\
			{bowtie2_map_parameters} \\
			--threads {NCORE} \\
			-x {TARGET_REF} \\
			-1 {INPUT1} \\
			-2 {INPUT2} \\
			--al-conc-gz {OUTPUT2} \\
			--un-conc-gz {OUTPUT3}  2>> {REPORT1} | \\
			$SAMTOOLS view \\
			--threads {NCORE} \\
			-Sh -f 2 /dev/stdin | \\
			$SAMTOOLS sort \\
			--threads {NCORE} \\
			-O bam -T {TEMP_PATH} \\
			-o {OUTPUT1} - \\
			>> {LOG_FILE} 2>&1
			
			$SAMTOOLS index \\
			-@ {NCORE} \\
			-b {OUTPUT1} \\
			>> {LOG_FILE} 2>&1

			$PICARD \\
			MarkDuplicates \\
			REMOVE_DUPLICATES=false \\
			I={OUTPUT1} \\
			O={OUTPUT1}.tmp \\
			M={REPORT2} \\
			>> {LOG_FILE} 2>&1

			$SAMTOOLS sort \\
			--threads {NCORE} \\
			-O bam -T {TEMP_PATH} \\
			-o {OUTPUT1}.nodup \\
			{OUTPUT1}.tmp \\
			>> {LOG_FILE} 2>&1

			$SAMTOOLS index \\
			-@ {NCORE} \\
			-b {OUTPUT1}.nodup \\
			>> {LOG_FILE} 2>&1

			$PICARD \\
			AddOrReplaceReadGroups \\
			RGID={SAMPLE} \\
			RGLB={SAMPLE} \\
			RGSM={SAMPLE} \\
			RGPL=illumina \\
			RGPU={SAMPLE} \\
			INPUT={OUTPUT1}.nodup \\
			OUTPUT={OUTPUT1}.nodup.tmp \\
			>> {LOG_FILE} 2>&1

			$SAMTOOLS sort \\
			--threads {NCORE} \\
			-O bam -T {TEMP_PATH} \\
			-o {OUTPUT1} \\
			{OUTPUT1}.nodup.tmp \\
			>> {LOG_FILE} 2>&1

			$SAMTOOLS index \\
			-@ {NCORE} \\
			-b {OUTPUT1} \\
			>> {LOG_FILE} 2>&1


			""".format(**IO_Dict)
		else:
			execution_script += """	
			$BOWTIE2 \\
			{bowtie2_map_parameters} \\
			--threads {NCORE} \\
			-x {TARGET_REF} \\
			-U {INPUT1} \\
			--al-gz {OUTPUT2} \\
			--un-gz {OUTPUT3}  2>> {REPORT1} | \\
			$SAMTOOLS view --threads {NCORE} -Sh  /dev/stdin | \\
			$SAMTOOLS sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1} - >> {LOG_FILE} 2>&1

			$SAMTOOLS index -@ {NCORE} -b {OUTPUT1} >> {LOG_FILE} 2>&1

			$PICARD \\
			MarkDuplicates \\
			REMOVE_DUPLICATES=false \\
			I={OUTPUT1} O={OUTPUT1}.tmp M={REPORT2} >> {LOG_FILE} 2>&1

			$SAMTOOLS sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1}.nodup {OUTPUT1}.tmp >> {LOG_FILE} 2>&1
			$SAMTOOLS index -@ {NCORE} -b {OUTPUT1}.nodup >> {LOG_FILE} 2>&1
			
			$PICARD \\
			AddOrReplaceReadGroups \\
			RGID={SAMPLE} \\
			RGLB={SAMPLE} \\
			RGSM={SAMPLE} \\
			RGPL=illumina \\
			RGPU={SAMPLE} \\
			INPUT={OUTPUT1}.nodup \\
			OUTPUT={OUTPUT1}.nodup.tmp \\
			>> {LOG_FILE} 2>&1

			$SAMTOOLS sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1} {OUTPUT1}.nodup.tmp >> {LOG_FILE} 2>&1
			$SAMTOOLS index -@ {NCORE} -b {OUTPUT1} >> {LOG_FILE} 2>&1
			
			""".format(**IO_Dict)
		execution_script += "# -------------------------------------------\n"

	else:
		##
		execution_script += extra_script
	
	return execution_script


def build_bowtie2_map_execution_script(general_Dict):
	"""
	"""
	execution_script = bowtie2_map_script(general_Dict)
	return execution_script
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #