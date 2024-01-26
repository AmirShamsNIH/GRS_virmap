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


def build_gatk_IO_dict(general_Dict):
	"""
	"""
	tag_List = [""]
	
	IO_Dict = {}
	# ++++++++++++++++++++++++++
	SNAKEFILE = general_Dict['SNAKEFILE']
	SNAKERULE = general_Dict['SNAKERULE']
	TITLE = general_Dict['TITLE']
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
	IO_Dict["TITLE"] = TITLE
	IO_Dict["DESIGN"] = DESIGN
	IO_Dict["INPUT_LAYOUT"] = general_Dict["METADATA"]["INPUT_LAYOUT"][0]
	IO_Dict["MAIN_PATH"] = general_Dict["MAIN_PATH"]
	IO_Dict["OUTPUT_PATH"] = general_Dict["OUTPUT_PATH"]
	IO_Dict["REPORT_PATH"] = general_Dict["REPORT_PATH"]
	IO_Dict["REFERENCE"] = general_Dict["METADATA"]["REFERENCE"]
	IO_Dict["TARGET_REF"] = general_Dict["METADATA"]["TARGET"]
	IO_Dict["CONTAMINATION"] = general_Dict["METADATA"]["CONTAMINATION"]
	# +++++	+++++++++++++++++++++
	#report
	IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
	IO_Dict["TEMP_PATH"] = general_Dict["TEMP_PATH"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++
	return IO_Dict


def build_gatk_haplotypecaller_script(IO_Dict):
	"""
	"""
	execution_script = """
		module load GATK
		module load samtools
		module load bcftools
		module load picard
		module load freebayes
	"""

	target_ref_path = IO_Dict["REFERENCE"]["TARGET_REF"]
	
	for each_input in IO_Dict["INPUT"]:
		"""
		"""
		IO_Dict["INPUT1"] = each_input
		bam_prefix = os.path.basename(each_input).replace(".bam", "")
		target = bam_prefix.split("bowtie2_map.")[1]
		IO_Dict["target_ref_file"] = target_ref_path + target + "/" + target + ".fa"

		IO_Dict["OUTPUT1"] = "{OUTPUT_PATH}".format(**IO_Dict) + bam_prefix + ".gatk_haplotypecaller.vcf.gz"
		IO_Dict["OUTPUT2"] = "{OUTPUT_PATH}".format(**IO_Dict) + bam_prefix + ".gatk_haplotypecaller.filtered.vcf.gz"

		IO_Dict["REPORT1"] = "{REPORT_PATH}".format(**IO_Dict) + bam_prefix + ".picard_nodup.txt"

		execution_script +="""
			freebayes -p {NCORE} --pooled-continuous --min-alternate-fraction 0.01 {INPUT1} -f {target_ref_file} \\
			> {OUTPUT1}.tmp

			cp {OUTPUT1}.tmp {OUTPUT2}.tmp 


			# gatk HaplotypeCaller \\
			# -ploidy 2 \\
			# -R {target_ref_file} -I {INPUT1} -O {OUTPUT1}.tmp >> {LOG_FILE} 2>&1
			# cp {OUTPUT1}.tmp {OUTPUT2}.tmp 
			# #bcftools filter \\
			# #-i "(TYPE="snp" && QUAL>500 && INFO/DP>20) || (TYPE="indel" && QUAL>1000 && INFO/DP>20)" -O v -o {OUTPUT2}.tmp {OUTPUT1}.tmp >> {LOG_FILE} 2>&1

			bgzip -c {OUTPUT1}.tmp > {OUTPUT1}
			bgzip -c {OUTPUT2}.tmp > {OUTPUT2}

			bcftools index {OUTPUT1}
			bcftools index {OUTPUT2}
			# cd {OUTPUT_PATH}
			
			# rm {OUTPUT1}.tmp*
			# rm {OUTPUT2}.tmp*
		
		""".format(**IO_Dict)

	return execution_script



def build_gatk_haplotypecaller_execution_script(general_Dict):
	"""
	"""
	
	IO_Dict = build_gatk_IO_dict(general_Dict)

	execution_script_String = build_gatk_haplotypecaller_script(IO_Dict)

	
	
	return execution_script_String

# ################################### CONFIGURATION ############################### #



def build_gatk_combine_script(IO_Dict):
	"""
	"""
	execution_script = """
		module load GATK
		module load samtools
		module load bcftools
	"""

	target_ref_path = IO_Dict["REFERENCE"]["TARGET_REF"]

	TARGET_REF = IO_Dict["TARGET_REF"]
	TARGTE_REF_STRING = "-".join(TARGET_REF)
	IO_Dict["TARGTE_REF_STRING"] = TARGTE_REF_STRING
	CONTAMINATION_REF = IO_Dict["CONTAMINATION"]
	CONTAMINATION_REF_STRING = "-".join(CONTAMINATION_REF)
	IO_Dict["CONTAMINATION_REF_STRING"] = CONTAMINATION_REF_STRING
	for each_target_ref in TARGET_REF:
		##
		filtered_input_List = []
		unfiltered_input_List = []

		IO_Dict["target_ref_file"] = target_ref_path + each_target_ref + "/" + each_target_ref + ".fa"
		
		IO_Dict["OUTPUT1"] = "{OUTPUT_PATH}{TITLE}.{DESIGN}.".format(**IO_Dict) + each_target_ref + ".aggregate.vcf.gz"

		IO_Dict["OUTPUT2"] = "{OUTPUT_PATH}{TITLE}.{DESIGN}.".format(**IO_Dict) + each_target_ref + ".aggregate.txt"
		

		for each_input in IO_Dict["INPUT"]:
			##
			
			filename = os.path.basename(each_input)
			if each_target_ref in filename:
				#
				if ".filtered.vcf.gz" in filename:
					filtered_input_List.append(each_input)
				else:
					unfiltered_input_List.append(each_input)
			else:
				pass
		else:
			pass
		#++++++++++++++++++++++++++
		execution_script +="""
			bcftools merge """ + " ".join(unfiltered_input_List) + """ > {OUTPUT1}.tmp 2>> {LOG_FILE}

			bgzip -c {OUTPUT1}.tmp > {OUTPUT1}
			
			rm -rf {OUTPUT1}.tmp
			
			gatk IndexFeatureFile -I {OUTPUT1} >> {LOG_FILE} 2>&1

			gatk VariantsToTable \\
			-F CHROM -F POS -F TYPE -F REF -F ALT -GF AD \\
			-R {target_ref_file} -V {OUTPUT1} -O {OUTPUT2} >> {LOG_FILE} 2>&1



			""".format(**IO_Dict)
	else:
		execution_script +="""
		python /data/shamsaddinisha/Development/GRS_virmap/library/collect_report.py {MAIN_PATH}../../ {OUTPUT_PATH} {TARGTE_REF_STRING} {CONTAMINATION_REF_STRING}
		""".format(**IO_Dict)

	return execution_script



def build_gatk_combine_execution_script(general_Dict):
	"""
	"""
	
	IO_Dict = build_gatk_IO_dict(general_Dict)

	execution_script_String = build_gatk_combine_script(IO_Dict)

	
	
	return execution_script_String



# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE RULES ############################## #
# ################################### FINITO ###################################### #