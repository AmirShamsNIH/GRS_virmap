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


def log_error(frameinfo):
	"""
	"""
	print("Error Occured at:", frameinfo.filename, frameinfo.lineno)
	sys.exit(2)
	return True

def build_cat_concat_execution_script(general_Dict):
	"""
	"""
	execution_script = ""
	the_index = 1
	for each_input in general_Dict["INPUT_DICT"]:
		##
		execution_script += """cat """ + general_Dict["INPUT" + the_index ] + """ >> """ + general_Dict["OUTPUT" + the_index ]
	else:
		##
		pass
	return [execution_script]


def build_execution_script(general_Dict):
	"""
	"""
	if general_Dict["INPUT_LAYOUT"] == "paired":
		bbtools_reformat_execution_script_List = build_paired_end_bbtools_reformat_execution_script(general_Dict)
	else:
		bbtools_reformat_execution_script_List = build_single_end_bbtools_reformat_execution_script(general_Dict)
	
	return [bbtools_reformat_execution_script_List]

# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #