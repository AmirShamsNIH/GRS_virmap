# ################################### INFO ######################################## #
# KRAKEN: trim fastq/fasta
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


def kraken2_map_script(general_Dict):
	"""
	https://github.com/DerrickWood/kraken2/blob/master/scripts/kraken2
	"""
	IO_Dict = {}
	execution_script = ""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	KRAKEN2_REFERENCE_LIST = general_Dict["CONFIG"]["REFERENCE"]["KRAKEN2_REFERENCE_INDEX"]
	PLATFORM = general_Dict["CONFIG"]["EXECUTION_PLATFORM"]
	REFERENECE_DICT = general_Dict["CONFIG"]["REFERENCE"]["KRAKEN2_REFERENCE"][PLATFORM]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if PLATFORM == "BIOWULF":
		execution_script += "module load kraken/2.1.2\n"
		execution_script += "module load kronatools\n"
		execution_script += "KRAKEN2=$(which kraken2)\n"
		execution_script += "KRONATOOLS=$(which ktImportTaxonomy)\n"
	elif PLATFORM == "BIGSKY":
		execution_script += "KRAKEN2=$(which kraken2)\n"
		execution_script += "KRONATOOLS=$(which ktImportTaxonomy)\n"
	execution_script += "# -------------------------------------------\n"

	for each_referenece in KRAKEN2_REFERENCE_LIST:
		##
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IO_Dict["KRAKEN2_REFERENCE"] = REFERENECE_DICT[each_referenece]
		# ----------------------------------------------------------------
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if len(general_Dict["INPUT"]) > 1:
			IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
			IO_Dict["INPUT2"] = general_Dict["INPUT"][1]
		else:
			IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
		#+
		IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_referenece + ".kraken2.txt"
		IO_Dict["REPORT2"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_referenece + ".kraken2_report.txt"
		IO_Dict["REPORT3"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_referenece + ".kraken2_report.html"
		#+
		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
		IO_Dict["SAMPLE"] = general_Dict["SAMPLE"]
		IO_Dict["kraken2_map_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["ALIGNMENT"]["kraken2_map"]
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
		if len(general_Dict["INPUT"]) > 1:
			#
			execution_script += """
				$KRAKEN2 \\
				--threads {NCORE} \\
				{kraken2_map_parameters} \\
				--db {KRAKEN2_REFERENCE} \\
				--report {REPORT2} \\
				{INPUT1} \\
				{INPUT2} \\
				--paired \\
				>> {REPORT1} \\
				2>> {LOG_FILE}
				
				$KRONATOOLS \\
				-q 2 -t 3 \\
				{REPORT2} \\
				-o {REPORT3} \\
				>> {LOG_FILE} 2>&1
				""".format(**IO_Dict)
		else:
			execution_script += """
				$KRAKEN2 \\
				--threads {NCORE} \\
				{kraken2_map_parameters} \\
				--db {KRAKEN2_REFERENCE} \\
				--report {REPORT2} \\
				{INPUT1} \\
				> {REPORT1} \\
				2>> {LOG_FILE}

				$KRONATOOLS \\
				-q 2 -t 3 \\
				{REPORT2} \\
				-o {REPORT3} \\
				>> {LOG_FILE} 2>&1
				""".format(**IO_Dict)
		
		execution_script += "# -------------------------------------------\n"

	else:
		##
		pass
	
	return execution_script


def build_kraken2_map_execution_script(general_Dict):
	"""
	"""
	execution_script = kraken2_map_script(general_Dict)
	return execution_script
# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #