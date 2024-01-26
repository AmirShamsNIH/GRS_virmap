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

def bbtools_reformat_script(general_Dict):
	"""
	https://github.com/BioInfoTools/BBMap/blob/master/sh/reformat.sh
	"""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	IO_Dict = {}
	
	IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
	IO_Dict["INPUT2"] = general_Dict["INPUT"][1]
	#+
	IO_Dict["OUTPUT1"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + ".R1.fastq.gz"
	IO_Dict["OUTPUT2"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + ".R2.fastq.gz"
	#+
	IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + ".txt"
	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	#+
	IO_Dict["NCORE"] = 5
	IO_Dict["bbtools_reformat_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["PRE_PROCESS"]["bbtools_reformat"]
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script = "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIOWULF":
		execution_script += "module load bbtools\n"
		execution_script += "BBTOOLS=\"$(which bbtools) reformat\"\n"
	elif general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIGSKY":
		execution_script += "BBTOOLS=$(which reformat.sh)\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# -------------------------------------------\n"
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"



	if len(general_Dict["INPUT"]) > 1:
		execution_script += """
		$BBTOOLS \\
		{bbtools_reformat_parameters} \\
		in={INPUT1} \\
		in2={INPUT2} \\
		out={OUTPUT1} \\
		out2={OUTPUT2} \\
		>> {LOG_FILE} 2>&1
		""".format(**IO_Dict)
	else:
		execution_script += """
		$BBTOOLS \\
		{bbtools_reformat_parameters} \\
		in={INPUT1} \\
		out={OUTPUT1} \\
		>> {LOG_FILE} 2>&1
		""".format(**IO_Dict)
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# -------------------------------------------\n"
	return execution_script


def build_bbtools_reformat_execution_script(general_Dict):
	"""
	"""

	execution_script = bbtools_reformat_script(general_Dict)

	return execution_script
# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #