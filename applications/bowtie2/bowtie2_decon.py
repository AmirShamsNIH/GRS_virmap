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
# ################################### BOWTIE2 DECON ############################### #

def bowtie2_decon_script(general_Dict):
	"""
	https://gensoft.pasteur.fr/docs/bowtie2/2.3.5.1/
	"""
	IO_Dict = {}
	execution_script = ""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	CONTAMINATION_LIST = general_Dict["CONFIG"]["REFERENCE"]["CONTAMINATION_REFERENCE_INDEX"]
	PLATFORM = general_Dict["CONFIG"]["EXECUTION_PLATFORM"]
	PREFIX = os.path.basename(general_Dict["INPUT"][0]).split(".R1.fastq.gz")[0]
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if PLATFORM == "BIOWULF":
		execution_script += "module load samtools\n"
		execution_script += "module load bowtie/2-2.5.1\n"
		execution_script += "SAMTOOLS=$(which samtools)\n"
		execution_script += "BOWTIE2=$(which bowtie2)\n"
	elif PLATFORM == "BIGSKY":
		execution_script += "SAMTOOLS=$(which samtools)\n"
		execution_script += "BOWTIE2=$(which bowtie2)\n"
	execution_script += "# -------------------------------------------\n"

	for each_contamination in CONTAMINATION_LIST:
		##
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		each_contamination_index  = CONTAMINATION_LIST.index(each_contamination)
		IO_Dict["CONTAMINATION_REF"] = general_Dict["CONFIG"]["REFERENCE"]["CONTAMINATION_REFERENCE"][PLATFORM][each_contamination]
		extra_script = "# +++++++++++++++++++++++++++++++++++++++++++\n"
		extra_script += "cp " + general_Dict["OUTPUT_PATH"] + PREFIX + "." + each_contamination + "_unmapped.R1.fastq.gz " + general_Dict["OUTPUT_PATH"] + PREFIX + ".bowtie2_decon.R1.fastq.gz \n"
		extra_script += "cp " + general_Dict["OUTPUT_PATH"] + PREFIX + "." + each_contamination + "_unmapped.R2.fastq.gz " + general_Dict["OUTPUT_PATH"] + PREFIX + ".bowtie2_decon.R2.fastq.gz \n"
		extra_script += "#rm -rf " + general_Dict["OUTPUT_PATH"] + PREFIX + "*_unmapped.R*.fastq.gz \n"
		extra_script += "# -------------------------------------------\n"
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if each_contamination_index == 0:
			#
			IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
			IO_Dict["INPUT2"] = general_Dict["INPUT"][1]
		else:
			#
			last_contamination_index = each_contamination_index - 1
			last_contamination = CONTAMINATION_LIST[last_contamination_index]
			IO_Dict["INPUT1"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + last_contamination + """_unmapped.R1.fastq.gz"""
			IO_Dict["INPUT2"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + last_contamination + """_unmapped.R2.fastq.gz"""
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IO_Dict["OUTPUT1"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + each_contamination + ".bam"
		IO_Dict["OUTPUT2"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + each_contamination + ".R%.fastq.gz"
		IO_Dict["OUTPUT3"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + each_contamination + "_unmapped.R%.fastq.gz"
		#+
		IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_contamination + ".txt"
		#+
		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
		IO_Dict["NCORE"] = 5
		IO_Dict["bowtie2_decon_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["PRE_PROCESS"]["bowtie2_decon"] 
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
		if len(general_Dict["INPUT"]) > 1:
			execution_script += """
			$BOWTIE2 \\
			--threads {NCORE} \\
			{bowtie2_decon_parameters} \\
			-x {CONTAMINATION_REF} \\
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
			""".format(**IO_Dict)
		else:
			execution_script += """
			$BOWTIE2 \\
			--threads {NCORE} \\
			{bowtie2_decon_parameters} \\
			-x {CONTAMINATION_REF} \\
			-U {INPUT[0]} \\
			--al-gz {OUTPUT2} \\
			--un-gz {OUTPUT3}  2>> {REPORT1} | \\
			$SAMTOOLS view --threads {NCORE} -Sh -F 4 /dev/stdin 2>> {LOG_FILE} | \\
			$SAMTOOLS sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1} - >> {LOG_FILE} 2>&1
			$SAMTOOLS index -@ {NCORE} -b {OUTPUT1} >> {LOG_FILE} 2>&1
			""".format(**IO_Dict)
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		execution_script += "# -------------------------------------------\n"

	else:
		##
		execution_script += extra_script
	
	return execution_script


def build_bowtie2_decon_execution_script(general_Dict):
	"""
	"""
	execution_script = bowtie2_decon_script(general_Dict)
	return execution_script

# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #