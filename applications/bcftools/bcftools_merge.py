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
def bcftools_merge_script(general_Dict):
	"""
	https://samtools.github.io/bcftools/bcftools.html
	"""
	IO_Dict = {}
	execution_script = ""
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	PLATFORM = general_Dict["CONFIG"]["EXECUTION_PLATFORM"]
	TARGET_REFERENCE_DICT = general_Dict["CONFIG"]["TARGET_REFERENCE"][PLATFORM]
	TARGET_REF = general_Dict["CONFIG"]["TARGET_REFERENCE_INDEX"]
	IO_Dict["SCRIPTS_PATH"] = general_Dict["CONFIG"]["WORKFLOW_BASEDIR"] + "/../scripts"
	
	TARGET_REF_STRING = "-".join(TARGET_REF)
	IO_Dict["TARGET_REF_STRING"] = TARGET_REF_STRING
	
	CONTAMINATION_REF = general_Dict["CONFIG"]["CONTAMINATION_REFERENCE_INDEX"]
	CONTAMINATION_REF_STRING = "-".join(CONTAMINATION_REF)
	IO_Dict["CONTAMINATION_REF_STRING"] = CONTAMINATION_REF_STRING
	IO_Dict["OUTPUT_PATH"] = general_Dict["OUTPUT_PATH"]
	IO_Dict["MAIN_PATH"] = general_Dict["MAIN_PATH"]
	IO_Dict["NCORE"] = 5
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	execution_script += "# +++++++++++++++++++++++++++++++++++++++++++\n"
	if PLATFORM == "BIOWULF":
		execution_script += "module load bcftools\n"
		execution_script += "module load freebayes\n"
		execution_script += "module load GATK\n"
		execution_script += "module load bcftools\n"
		execution_script += "module load rtg-tools\n"
		execution_script += "module load snpEff\n"
		execution_script += "FREEBAYES=$(which freebayes)\n"
		execution_script += "BCFTOOLS=$(which bcftools)\n"
		execution_script += "GATK=$(which gatk)\n"
		execution_script += "BGZIP=$(which bgzip)\n"
		execution_script += "RTG=$(which rtg)\n"
		execution_script += "SNPEFF=\"java -Xmx12g -jar $SNPEFF_JAR\"\n"
	elif PLATFORM == "BIGSKY":
		execution_script += "FREEBAYES=$(which freebayes)\n"
		execution_script += "BCFTOOLS=$(which bcftools)\n"
		execution_script += "GATK=$(which gatk)\n"
		execution_script += "BGZIP=$(which bgzip)\n"
		execution_script += "RTG=/gs1/RTS/NextGen/bin/grs_virmap/bin/rtg\n"
		execution_script += "SNPEFF=\"java -jar /gs1/RTS/NextGen/snpEff/snpEff.jar\"\n"
	execution_script += "# -------------------------------------------\n"
	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for each_target_ref in TARGET_REFERENCE_DICT:
		##
		vcf_List = []
		IO_Dict["OUTPUT1"] = general_Dict["OUTPUT_PATH"]  + general_Dict["TITLE"] + "." + general_Dict["DESIGN"] + "." + each_target_ref + ".aggregate.vcf.gz"
		IO_Dict["OUTPUT2"] = general_Dict["OUTPUT_PATH"]  + general_Dict["TITLE"] + "." + general_Dict["DESIGN"] + "." + each_target_ref + ".aggregate.txt"
		IO_Dict["OUTPUT3"] = general_Dict["OUTPUT_PATH"]  + general_Dict["TITLE"] + "." + general_Dict["DESIGN"] + "." + each_target_ref + ".aggregate.snpeff.vcf.gz"
		#+
		IO_Dict["REPORT1"] = general_Dict["REPORT_PATH"] + general_Dict["TITLE"] + "." + general_Dict["DESIGN"] + "." + each_target_ref + ".aggregate.bcftools_stats.txt"
		IO_Dict["REPORT2"] = general_Dict["REPORT_PATH"] + general_Dict["TITLE"] + "." + general_Dict["DESIGN"] + "." + each_target_ref + ".aggregate.rtg_vcfstats.txt"
		IO_Dict["REPORT3"] = general_Dict["REPORT_PATH"] + general_Dict["TITLE"] + "." + general_Dict["DESIGN"] + "." + each_target_ref + ".aggregate.snpeff_summary.csv"
		IO_Dict["REPORT4"] = general_Dict["REPORT_PATH"] + general_Dict["TITLE"] + "." + general_Dict["DESIGN"] + "." + each_target_ref + ".aggregate.snpeff_summary.html"
		#+
		IO_Dict["TARGET_REF_NAME"] = each_target_ref
		#IO_Dict["TARGET_REF"] = TARGET_REFERENCE_DICT[each_target_ref]

		OUTPUT = general_Dict["CONFIG"]["OUTPUT"]
		TITLE = general_Dict["CONFIG"]["TITLE"]
		IO_Dict["TARGET_REF"] = OUTPUT + "/" + TITLE + "/custom_virmapDB/output/" + each_target_ref + "/" + each_target_ref + ".fa"
		IO_Dict["SNPEFF_CONFIG"] = OUTPUT + "/" + TITLE + "/custom_virmapDB/output/snpEff.config"

		IO_Dict["LOG_FILE"] = general_Dict["LOG_FILE"]
		IO_Dict["NCORE"] = 5

		IO_Dict["bcftools_merge_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["REPORT"]["bcftools_merge"]
		IO_Dict["snpeff_ann_parameters"] = general_Dict["CONFIG"]["PARAMETERS"]["REPORT"]["snpeff_ann"]
		for each_input in general_Dict["INPUT"]:
			##
			filename = os.path.basename(each_input)
			if each_target_ref in filename:
				#
				if "snpeff" in filename:
					continue

				vcf_List.append(each_input)
			else:
				#
				pass
			
		else:
			##
			
			#++++++++++++++++++++++++++
			execution_script +="""
				$BCFTOOLS merge  \\
				""" + " ".join(vcf_List) + """ > {OUTPUT1}.tmp >> \\
				{LOG_FILE} 2>&1

				$BGZIP -c {OUTPUT1}.tmp > {OUTPUT1}
				
				rm -rf {OUTPUT1}.tmp


				$BCFTOOLS index -f  \\
				{OUTPUT1} >> \\
				{LOG_FILE} 2>&1
				
				$GATK IndexFeatureFile \\
				-I {OUTPUT1} \\
				>> {LOG_FILE} 2>&1

				$GATK VariantsToTable \\
				-F CHROM -F POS -F TYPE -F REF -F ALT -GF AD \\
				-R {TARGET_REF} \\
				-V {OUTPUT1} \\
				-O {OUTPUT2} \\
				>> {LOG_FILE} 2>&1

				$BCFTOOLS stats {OUTPUT1} \\
				> {REPORT1} 2>> {LOG_FILE}

				$RTG vcfstats {OUTPUT1} \\
				> {REPORT2} 2>> {LOG_FILE}

				$SNPEFF ann \\
				{snpeff_ann_parameters} \\
				-config {SNPEFF_CONFIG} -v {TARGET_REF_NAME} -csvStats {REPORT3} -stats {REPORT4} \\
				{OUTPUT1} > {OUTPUT3} 2>> {LOG_FILE} 

			""".format(**IO_Dict)
			#++++++++++++++++++++++++++
	else:
		execution_script +="""
		python {SCRIPTS_PATH}/collect_report.py {MAIN_PATH}../../ {OUTPUT_PATH} {TARGET_REF_STRING} {CONTAMINATION_REF_STRING}
		""".format(**IO_Dict)

	return execution_script
   

def build_bcftools_merge_execution_script(general_Dict):
	"""
	"""
	execution_script = bcftools_merge_script(general_Dict)
	return execution_script