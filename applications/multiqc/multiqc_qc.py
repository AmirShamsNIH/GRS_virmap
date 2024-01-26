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

def multiqc_qc_script(general_Dict):
	"""
	https://multiqc.info/docs/getting_started/running_multiqc/
	"""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	IO_Dict = {}
	IO_Dict["SAMPLE"] = general_Dict["SAMPLE"]
	IO_Dict["TITLE"] = general_Dict["TITLE"]
	IO_Dict["DESIGN"] = general_Dict["DESIGN"]
	IO_Dict["MAIN_PATH"] = general_Dict["MAIN_PATH"]
	IO_Dict["REPORT_PATH"] = general_Dict["REPORT_PATH"]
	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script = ""
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIOWULF":
		execution_script += "module load multiqc\n"
		execution_script += "MULTIQC=$(which multiqc)\n"
	elif general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIGSKY":
		execution_script += "MULTIQC=$(which multiqc)\n"
	execution_script += "# -------------------------------------------\n"
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	execution_script += """
		$MULTIQC \\
		--force --exclude general_stats \\
		--no-ansi \\
		-s \\
		--outdir {REPORT_PATH} \\
		--filename {TITLE}.{DESIGN}.{SAMPLE} \\
		{MAIN_PATH} \\
		>> {LOG_FILE} 2>&1
	""".format(**IO_Dict)
	execution_script += "# -------------------------------------------\n"

	return execution_script


def build_multiqc_qc_execution_script(general_Dict):
	"""
	"""

	execution_script = multiqc_qc_script(general_Dict)

	return execution_script
# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #