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

def fastp_trim_script(general_Dict):
	"""
	https://github.com/OpenGene/fastp#input-and-output
	"""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	IO_Dict = {}

	IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
	IO_Dict["INPUT2"] = general_Dict["INPUT"][1]
	PREFIX = os.path.basename(general_Dict["INPUT"][0]).split(".R1.fastq.gz")[0]
	#+
	IO_Dict["OUTPUT1"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + general_Dict["SNAKERULE"] + ".R1.fastq.gz"
	IO_Dict["OUTPUT2"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + general_Dict["SNAKERULE"] + ".R2.fastq.gz"
	IO_Dict["OUTPUT3"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + "fastp_discarded" + ".R1.fastq.gz"
	IO_Dict["OUTPUT4"] = general_Dict["OUTPUT_PATH"] + PREFIX + "." + "fastp_discarded" + ".R2.fastq.gz"
	#+
	IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + ".fastp.html"
	IO_Dict["REPORT2"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + ".fastp.json"
	#+
	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	IO_Dict["NCORE"] = 5
	IO_Dict["fastp_trim_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["PRE_PROCESS"]["fastp_trim"] 
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script = "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIOWULF":
		execution_script += "module load fastp\n"
		execution_script += "FASTP=$(which fastp)\n"
	elif general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIGSKY":
		execution_script += "FASTP=$(which /gs1/RTS/NextGen/bin/grs_virmap/bin/fastp)\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# -------------------------------------------\n"
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if len(general_Dict["INPUT"]) > 1:
		execution_script += """
		$FASTP \\
		--thread {NCORE} \\
		{fastp_trim_parameters} \\
		--detect_adapter_for_pe \\
		--in1 {INPUT1} \\
		--in2 {INPUT2} \\
		--out1 {OUTPUT1} \\
		--out2 {OUTPUT2} \\
		--unpaired1 {OUTPUT3} \\
		--unpaired2 {OUTPUT4} \\
		--html {REPORT1} \\
		--json {REPORT2} \\
		>> {LOG_FILE} 2>&1
		""".format(**IO_Dict)
	else:
		execution_script += """
		$FASTP \\
		--thread {NCORE} \\
		{fastp_trim_parameters} \\
		--in1 {INPUT1} \\
		--out1 {OUTPUT1} \\
		--failed_out {OUTPUT3} \\
		--html {REPORT1} \\
		--json {REPORT2} \\
		>> {LOG_FILE} 2>&1
		""".format(**IO_Dict)
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# -------------------------------------------\n"
	return execution_script


def build_fastp_trim_execution_script(general_Dict):
	"""
	"""
	execution_script = fastp_trim_script(general_Dict)
	return execution_script

# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #