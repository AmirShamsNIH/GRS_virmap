# ################################### INFO ######################################## #
# snpeff: trim fastq/fasta
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


def build_snpeff_annotate_IO_dict(general_Dict):
	"""
	"""
	
	IO_Dict = {}
	# ++++++++++++++++++++++++++
	SNAKEFILE = general_Dict['SNAKEFILE']
	SNAKERULE = general_Dict['SNAKERULE']
	DESIGN = general_Dict['DESIGN']
	SAMPLE = general_Dict['SAMPLE']
	# ++++++++++++++++++++++++++
	#Input
	item_List = []
	snakerule_input_List = general_Dict["IO"][SNAKEFILE][SNAKERULE][DESIGN][SAMPLE]["input"]
	for each_input in snakerule_input_List:
		##
		item_List.extend(glob.glob(each_input))
	else:
		pass
	item_List = list(set(item_List))
	item_List.sort()
	IO_Dict["INPUT"] = item_List
	# ++++++++++++++++++++++++++
	# ++++++++++++++++++++++++++
	IO_Dict["SAMPLE"] = SAMPLE
	IO_Dict["INPUT_LAYOUT"] = general_Dict["METADATA"]["INPUT_LAYOUT"][0]
	IO_Dict["OUTPUT_PATH"] = general_Dict["OUTPUT_PATH"]
	IO_Dict["REPORT_PATH"] = general_Dict["REPORT_PATH"] 
	IO_Dict["TARGET"] = general_Dict["METADATA"]["TARGET"]
	IO_Dict["contamination"] = general_Dict["METADATA"]["CONTAMINATION"]
	IO_Dict["REFERENCE"] = general_Dict["METADATA"]["REFERENCE"]
	# +++++	+++++++++++++++++++++
	#report
	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++
	return IO_Dict


def build_kraken2_viral_script(IO_Dict):
	"""
	"""
	execution_script = """
		module load kraken/2.1.2
		module load kronatools
	"""
	
	

	IO_Dict["TARET_REF_FILE"] = IO_Dict["REFERENCE"]["KRAKEN2_VIRAL"]
	IO_Dict["INPUT1"] = "{INPUT[0]}".format(**IO_Dict)
	IO_Dict["INPUT2"] = "{INPUT[1]}".format(**IO_Dict)	
	#-
	#+
	# IO_Dict["OUTPUT1"] = "{OUTPUT_PATH}{SAMPLE}".format(**IO_Dict) + ".kraken2_viral." + each_target + ".bam"
	# IO_Dict["OUTPUT2"] = "{OUTPUT_PATH}{SAMPLE}".format(**IO_Dict) + ".kraken2_viral." + each_target + "_mapped.R%.fastq.gz"
	# IO_Dict["OUTPUT3"] = "{OUTPUT_PATH}{SAMPLE}".format(**IO_Dict) + ".kraken2_viral." + each_target + "_unmapped.R%.fastq.gz"
	#-
	#+
	IO_Dict["REPORT1"] = "{REPORT_PATH}{SAMPLE}".format(**IO_Dict) + ".kraken2_viral.txt"
	IO_Dict["REPORT2"] = "{REPORT_PATH}{SAMPLE}".format(**IO_Dict) + ".kraken2_viral.html"

	execution_script += """
	kraken2 \\
	--threads {NCORE} \\
	--fastq-input \\
	--db {TARET_REF_FILE} \\
	--report {REPORT1} --gzip-compressed \\
	{INPUT1} \\
	{INPUT2} \\
	--paired >> {LOG_FILE} 2>&1
	
	ktImportTaxonomy -q 2 -t 3 {REPORT1} \\
	-o {REPORT2} >> {LOG_FILE} 2>&1
	""".format(**IO_Dict)
	
	return execution_script	

def build_snpeff_annotate_script(IO_Dict):
	"""
	"""
	execution_script = """
		module load snpeff
	"""

	execution_script += """
	java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR \\
	-no-intergenic -no-upstream -no-downstream -no-utr \\
	
	kraken2 \\
	--threads {NCORE} \\
	--fastq-input \\
	--db {TARET_REF_FILE} \\
	--report {REPORT1} --gzip-compressed \\
	{INPUT1} \\
	{INPUT2} \\
	--paired >> {LOG_FILE} 2>&1
	
	ktImportTaxonomy -q 2 -t 3 {REPORT1} \\
	-o {REPORT2} >> {LOG_FILE} 2>&1
	""".format(**IO_Dict)
	
	return execution_script	
	
	

	IO_Dict["TARET_REF_FILE"] = IO_Dict["REFERENCE"]["KRAKEN2_VIRAL"]
	IO_Dict["INPUT1"] = "{INPUT[0]}".format(**IO_Dict)
	IO_Dict["INPUT2"] = "{INPUT[1]}".format(**IO_Dict)	

def build_kraken2_viral_execution_script(general_Dict):
	"""
	"""
	IO_Dict = {}
	INPUT_LAYOUT = general_Dict["METADATA"]["INPUT_LAYOUT"][0]
	IO_Dict = build_kraken2_viral_IO_dict(general_Dict)
	
	execution_script_List = build_kraken2_viral_script(IO_Dict)
	
	return execution_script_List


def build_snpeff_annotate_execution_script(general_Dict):
	"""
	"""
	IO_Dict = build_snpeff_annotate_IO_dict(general_Dict)
	
	execution_script_List = build_snpeff_annotate_script(IO_Dict)

	return execution_script_List

# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #