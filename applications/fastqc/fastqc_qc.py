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


def fastqc_qc_script(general_Dict):
	"""
	https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help
	"""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	IO_Dict = {}

	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	IO_Dict["REPORT_PATH"] = general_Dict["REPORT_PATH"]
	IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script = ""
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIOWULF":
		execution_script += "module load fastqc/0.11.9\n"
		execution_script += "FASTQC=$(which fastqc)\n"
	elif general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIGSKY":
		execution_script += "FASTQC=$(which fastqc)\n"
	execution_script += "# -------------------------------------------\n"
	for each_input in general_Dict["INPUT"]:
		##
		IO_Dict["INPUT1"] = each_input
		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		IO_Dict["REPORT_PATH"] = general_Dict["REPORT_PATH"]
		IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
		IO_Dict["NCORE"] = 5
		execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		execution_script += """
		$FASTQC \\
		--outdir {REPORT_PATH} \\
		--threads {NCORE} \\
		--dir {TEMP_PATH} \\
		{INPUT1} \\
		>> {LOG_FILE} 2>&1
		""".format(**IO_Dict)
		execution_script += "# -------------------------------------------\n"

	else:
		##
		pass

	return execution_script


def build_fastqc_qc_execution_script(general_Dict):
	"""
	"""

	execution_script = fastqc_qc_script(general_Dict)

	return execution_script
# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #