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


def build_picard_IO_dict(general_Dict):
	"""
	"""
	tag_List = [""]
	
	IO_Dict = {}
	# ++++++++++++++++++++++++++
	SNAKEFILE = general_Dict['SNAKEFILE']
	SNAKERULE = general_Dict['SNAKERULE']
	DESIGN = general_Dict['DESIGN']
	SAMPLE = general_Dict['SAMPLE']
	# ++++++++++++++++++++++++++
	#Input
	item_List = []
	snakerule_input_List = general_Dict["METADATA"]["snakefile"][SNAKEFILE][SNAKERULE][DESIGN][SAMPLE]["input"]
	for each_input in snakerule_input_List:
		##
		item_List.extend(glob.glob(each_input))
	else:
		pass
	processed_item_List = []
	for each_tag in tag_List:
		for each_item in item_List:
			if each_tag in os.path.basename(each_item):
				processed_item_List.append(each_item)
			else:
				pass
	item_List = list(set(processed_item_List))
	item_List.sort()
	IO_Dict["INPUT"] = item_List
	# ++++++++++++++++++++++++++
	# ++++++++++++++++++++++++++
	IO_Dict["SAMPLE"] = SAMPLE
	IO_Dict["INPUT_LAYOUT"] = general_Dict["METADATA"]["INPUT_LAYOUT"][0]
	IO_Dict["OUTPUT_PATH"] = general_Dict["OUTPUT_PATH"]
	IO_Dict["REPORT_PATH"] = general_Dict["REPORT_PATH"]
	# +++++	+++++++++++++++++++++
	#report
	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++
	return IO_Dict


def build_picard_nodup_script(IO_Dict):
	"""
	"""
	execution_script = """
		module load picard
		module load samtools
	"""
	
	for each_input in IO_Dict["INPUT"]:
		"""
		"""
		
		IO_Dict["INPUT1"] = each_input
		bam_prefix = os.path.basename(each_input).replace(".bam", "")
		IO_Dict["OUTPUT1"] = "{OUTPUT_PATH}".format(**IO_Dict) + bam_prefix + ".picard_nodup.bam"
		IO_Dict["REPORT1"] = "{REPORT_PATH}".format(**IO_Dict) + bam_prefix + ".picard_nodup.txt"
		execution_script +="""
			java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\
			MarkDuplicates \\
			REMOVE_DUPLICATES=true \\
			I={INPUT1} O={OUTPUT1}.tmp M={REPORT1} >> {LOG_FILE} 2>&1

			samtools sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1}.nodup {OUTPUT1}.tmp >> {LOG_FILE} 2>&1
			samtools index -@ {NCORE} -b {OUTPUT1}.nodup >> {LOG_FILE} 2>&1

			java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar \\
			AddOrReplaceReadGroups \\
			RGID={SAMPLE} RGLB={SAMPLE} RGSM={SAMPLE} RGPL=illumina RGPU={SAMPLE} INPUT={OUTPUT1}.nodup OUTPUT={OUTPUT1}.tmp >> {LOG_FILE} 2>&1

			samtools sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1} {OUTPUT1}.tmp >> {LOG_FILE} 2>&1
			samtools index -@ {NCORE} -b {OUTPUT1} >> {LOG_FILE} 2>&1
		
		""".format(**IO_Dict)

	return execution_script



def build_picard_nodup_execution_script(general_Dict):
	"""
	"""
	
	IO_Dict = build_picard_IO_dict(general_Dict)

	execution_script_String = build_picard_nodup_script(IO_Dict)
	
	return execution_script_String




def picard_markduplicate_script(general_Dict):
	"""
	https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
	"""
	IO_Dict = {}
	execution_script = ""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	target_List = general_Dict["METADATA"]["TARGET"]
	target_ref_path = general_Dict["METADATA"]["REFERENCE"]["TARGET_REF"]
	
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# ++++++++++++++++++++++++++++++++++++++++++++\n"
	if general_Dict["METADATA"]["PLATFORM"][0] == "BIOWULF":
		execution_script += "module load samtools\n"
		execution_script += "module load bowtie/2-2.5.1\n"
		execution_script += "samtools=$(which samtools)\n"
		execution_script += "bowtie2=$(which bowtie2)\n"
	elif general_Dict["METADATA"]["PLATFORM"][0] == "BIGSKY":
		execution_script += "samtools=$(which samtools)\n"
		execution_script += "bowtie2=$(which bowtie2)\n"
	execution_script += "# ++++++++++++++++++++++++++++++++++++++++++++\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for each_target in target_List:
		##
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		target_ref_file = target_ref_path + each_target + "/" + each_target + ".fa"
		IO_Dict["TARET_REF_FILE"] = target_ref_file
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if len(general_Dict["INPUT"]) > 1:
			IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
			IO_Dict["INPUT2"] = general_Dict["INPUT"][1]
		else:
			IO_Dict["INPUT1"] = general_Dict["INPUT"][0]
		IO_Dict["OUTPUT1"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." +  each_target + ".bam"
		IO_Dict["OUTPUT2"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_target + "_mapped.R%.fastq.gz"
		IO_Dict["OUTPUT3"] = general_Dict["OUTPUT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_target + "_unmapped.R%.fastq.gz"
		IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + "." + each_target + ".txt"
		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
		IO_Dict["NCORE"] = 5
		# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if len(general_Dict["INPUT"]) > 1:
			execution_script += """	
			bowtie2 \\
			-k10 -X 1500 \\
			--threads {NCORE} \\
			-x {TARET_REF_FILE} \\
			-1 {INPUT1} \\
			-2 {INPUT2} \\
			--al-conc-gz {OUTPUT2} \\
			--un-conc-gz {OUTPUT3}  2>> {REPORT1} | \\
			samtools view --threads {NCORE} -Sh -f 2 -F 1840 /dev/stdin | \\
			samtools sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1} - >> {LOG_FILE} 2>&1

			samtools index -@ {NCORE} -b {OUTPUT1} >> {LOG_FILE} 2>&1
			""".format(**IO_Dict)
		else:
			execution_script += """	
			bowtie2 \\
			-k10 \\
			--threads {NCORE} \\
			-x {TARET_REF_FILE} \\
			-U {INPUT1} \\
			--al-gz {OUTPUT2} \\
			--un-gz {OUTPUT3}  2>> {REPORT1} | \\
			samtools view --threads {NCORE} -Sh -F 260 /dev/stdin | \\
			samtools sort --threads {NCORE} -O bam -T {TEMP_PATH} -o {OUTPUT1} - >> {LOG_FILE} 2>&1

			samtools index -@ {NCORE} -b {OUTPUT1} >> {LOG_FILE} 2>&1
			""".format(**IO_Dict)
		execution_script += "# -------------------------------------------\n"

	else:
		##
		pass
	
	return execution_script


def build_picard_markduplicate_execution_script(general_Dict):
	"""
	"""
	execution_script = picard_markduplicate_script(general_Dict)
	return [execution_script]

# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #