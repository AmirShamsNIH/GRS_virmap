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


def samtools_idxstats_script(general_Dict):
	"""
	https://www.htslib.org/doc/samtools.html
	"""
	IO_Dict = {}
	execution_script = ""
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIOWULF":
		execution_script += "module load samtools\n"
		execution_script += "SAMTOOLS=$(which samtools)\n"
	elif general_Dict["CONFIG"]["EXECUTION_PLATFORM"] == "BIGSKY":
		execution_script += "SAMTOOLS=$(which samtools)\n"
	execution_script += "# -------------------------------------------\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for each_input in general_Dict["INPUT"]:
		##
		IO_Dict["INPUT1"] = each_input
		PREFIX = os.path.basename(each_input).replace(".bam", "")
		IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + PREFIX + ".samtools_idxstats.txt"
		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		IO_Dict["NCORE"] = 5
		execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
		execution_script +="""
			$SAMTOOLS \\
			idxstats \\
			--threads {NCORE} \\
			{INPUT1} > {REPORT1} \\
			2> {LOG_FILE}
		""".format(**IO_Dict)
		execution_script += "# --------------------------------------------\n"
	else:
		##
		pass
		
	return execution_script


def build_samtools_idxstats_execution_script(general_Dict):
	"""
	"""
	execution_script = samtools_idxstats_script(general_Dict)
	return execution_script


# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #