# ################################### INFO ######################################## #
# Author: Amir Shams
# Project: a project
# Aim: ...
# ################################### IMPORT ###################################### #

import os
import sys
import re
import glob
import json
import pandas
import importlib
import shutil
from datetime import datetime
from inspect import currentframe, getframeinfo
# ################################### FUNCTIONS ################################### #
DEBUG_MODE = False

def log_error(frameinfo):
	"""
	"""
	print("Error Occured at:", frameinfo.filename, frameinfo.lineno)
	sys.exit(2)
	return True

def is_file_exist(file_Path):
	"""
	"""
	if os.path.isfile(file_Path) and os.path.exists(file_Path) and os.access(file_Path, os.R_OK):
		return True
	else:
		return False

def add_date_tag(file_path):
	"""
	"""
	if is_file_exist(file_path):
		#
		current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
		str_current_datetime = str(current_datetime)
		file_name = os.path.basename(file_path)
		file_suffix = file_name.split(".")[-1]
		new_file_name = file_name.replace(file_suffix, str_current_datetime + "." + file_suffix)
		new_file_path = file_path.replace(file_name, new_file_name)
		shutil.copyfile(file_path, new_file_path)
		open(file_path, 'w').close()
	else:
		#
		pass
	return True

def build_metadata_dict(metadata_file_Path):
	##
	#metadata_file_Path = config_Dict["metadata"]
	metadata_Dict = {}
	#
	metadata_DF = pandas.read_csv(
		metadata_file_Path,
		encoding=None,
		skip_blank_lines=True,
		delimiter=","
	)
	metadata_DF.columns = metadata_DF.columns.str.lower()
	#
	metadata_List = metadata_DF.fillna('empty').to_dict('records')
	metadata_Dict = {}
	for each_item in metadata_List:
		##
		design = each_item["design"]
		sample = each_item["sample"]
		if design not in metadata_Dict:
			#
			metadata_Dict[design] = {}
			metadata_Dict[design][sample] = {**each_item}
		else:
			metadata_Dict[design][sample] = {**each_item}

	else:
		##
		pass
	
	return metadata_Dict


def get_last_snakerule(config_Dict, snakefile, last_snakefile, snakerule, app_Dict):
	"""
	"""
	app_List = app_Dict[snakefile]
	if app_List[0] == snakerule:
		#this is the first rule
		if last_snakefile == None:
			last_snakerule = None
		else:
			last_app_List = app_Dict[last_snakefile]
			last_snakerule = last_app_List[-1]
	else:
		snakerule_index = app_List.index(snakerule)
		last_snakerule = app_List[snakerule_index - 1]
	
	return last_snakerule


def get_app_json_dict(config_Dict, snakerule):
	"""
	"""
	# ++++++++++++++++++++++++++++++
	workflow_basedir = config_Dict["WORKFLOW_BASEDIR"]
	application_folder = snakerule.split("_")[0] #require more unified method
	application_json_file = workflow_basedir + "/../applications/" + application_folder + "/" + application_folder + ".json"
	
	with open(application_json_file) as f:
		application_json_Dict = json.load(f)
	# ++++++++++++++++++++++++++++++
	if snakerule in application_json_Dict:
		#
		#app_json_Dict = application_json_Dict[snakerule]
		return application_json_Dict
	else:
		print("snakerule: " + snakerule + " is not available in app json")
		print("application_json_file: " + application_json_file)
		log_error(getframeinfo(currentframe()))
	return False


def process_paired_fastq(forward_fastq_List, reverse_fastq_List):
	"""
	"""
	fastq_extension_regex = "(\.fastq\.gz$|\.fq\.gz$|\.fastq$|\.fq$)"
	forward_name_tag_regex = "(_R1|\.r1|\.R1)"
	reverse_name_tag_regex = "(_R2|\.r2|\.R2)"
	
	forward_fastq_List = list(set(forward_fastq_List))
	reverse_fastq_List = list(set(reverse_fastq_List))
	forward_fastq_List.sort()
	reverse_fastq_List.sort()
	
	processed_forward_fastq_List = []
	processed_reverse_fastq_List = []
	for each_forward_fastq, each_reverse_fastq in zip(forward_fastq_List, reverse_fastq_List):
		##
		forward_fastq_name = os.path.basename(each_forward_fastq)
		forward_fastq_name_no_extension = re.sub(fastq_extension_regex, "", forward_fastq_name)
		forward_fastq_name_no_extension_no_tag = re.sub(forward_name_tag_regex, "", forward_fastq_name_no_extension)
		# ++++++++++
		reverse_fastq_name = os.path.basename(each_reverse_fastq)
		reverse_fastq_name_no_extension = re.sub(fastq_extension_regex, "", reverse_fastq_name)
		reverse_fastq_name_no_extension_no_tag = re.sub(reverse_name_tag_regex, "", reverse_fastq_name_no_extension)
		# ++++++++++
		
		if(forward_fastq_name_no_extension_no_tag == reverse_fastq_name_no_extension_no_tag):
			processed_forward_fastq_List.append(each_forward_fastq)
			processed_reverse_fastq_List.append(each_reverse_fastq)
		else:
			#
			pass
	else:
		##
		pass
	
	if len(processed_forward_fastq_List) == 0 :
		#
		print("no paired-end fastq detected")
		log_error(getframeinfo(currentframe()))
	else:
		#
		return(processed_forward_fastq_List, processed_reverse_fastq_List)


def get_fastq(sample, input_List, fastq_layout):
	"""
	"""
	fastq_extension_List = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
	fastq_extension_regex = "(\.fastq\.gz$|\.fq\.gz$|\.fastq$|\.fq$)"
	forward_name_tag_List = ["_R1", ".r1", ".R1"]
	reverse_name_tag_List = ["_R2", ".r2", ".R2"]
	if fastq_layout == "paired":
		#
		forward_fastq_List = []
		reverse_fastq_List = []
		for each_input in input_List:
			##
			for each_fastq_extension in fastq_extension_List:
				##
				# ++++++++++++++++++++++
				for each_forward_name_tag in forward_name_tag_List:
					##
					forward_fastq_List.extend(glob.glob(each_input + "/" + sample + "*" + each_forward_name_tag + "*" + each_fastq_extension, recursive=True))
				else:
					##for each_forward_name_tag in forward_name_tag_List:
					pass
				# ++++++++++++++++++++++
				for each_reverse_name_tag in reverse_name_tag_List:
					##
					reverse_fastq_List.extend(glob.glob(each_input + "/" + sample + "*" + each_reverse_name_tag + "*" + each_fastq_extension, recursive=True))
				else:
					##for each_reverse_name_tag in reverse_name_tag_List:
					pass
				# ++++++++++++++++++++++
			else:
				##for each_fastq_extension in fastq_extension_List:
				pass
		else:
			##for each_input in input_List:
			pass
		# ++++++++++++++++++++++
		
		processed_forward_fastq_List, processed_reverse_fastq_List = process_paired_fastq(forward_fastq_List, reverse_fastq_List)
		return (processed_forward_fastq_List, processed_reverse_fastq_List)

	elif fastq_layout == "single":
		#
		forward_fastq_List = []
		reverse_fastq_List = []
		for each_input in input_List:
			##
			for each_fastq_extension in fastq_extension_List:
				##
				forward_fastq_List.extend(glob.glob(each_input + "/" + sample + "*" +  each_fastq_extension, recursive=True))
				reverse_fastq_List.append("")
			else:
				##for each_fastq_extension in fastq_extension_List:
				pass
		else:
			##for each_input in input_List:
			pass
		if len(forward_fastq_List)==0:
			#
			print("no fastq detected!")
			log_error(getframeinfo(currentframe()))
		else:

			return (forward_fastq_List, reverse_fastq_List)
	else:
		#
		print("OOPS Error!")
		log_error(getframeinfo(currentframe()))
		sys.exit(2)
	return False

# ################################### SNAKEMAKE ###########################################


def build_snakerule_input_list(config_Dict, design, sample, snakefile, snakerule, last_snakerule_List):
	"""
	"""
	# ++++++++++++++++++++++++++++++
	TITLE = config_Dict["TITLE"]
	INPUT_LIST = config_Dict["INPUT"]
	OUTPUT = config_Dict["OUTPUT"]
	# ++++++++++++++++++++++++++++++
	snakerule_json_Dict = get_app_json_dict(config_Dict, snakerule)
	# ++++++++++++++++++++++++++++++++
	snakerule_input_List = []
	for last_snakerule in last_snakerule_List:
		##
		if last_snakerule == "start_point":
			#this is starting point
			for each_input_extension in snakerule_json_Dict[snakerule]["input_extension"]:
				##
				for each_input in INPUT_LIST:
					##
					snakerule_input_List.append(each_input + "/" + sample + "*" + each_input_extension)
				else:
					##
					pass
			else:
				##
				pass
		elif last_snakerule == "None":
			#
			snakerule_input_List = []

		elif last_snakerule == "*":
			#this is snakefile sample report
			for each_input_extension in snakerule_json_Dict[snakerule]["input_extension"]:
				snakerule_input_List.append(OUTPUT + TITLE + "/" + design + "/" + snakefile + "/" + sample + "/**/" + sample + "*" +  each_input_extension)
			else:
				pass
		elif last_snakerule.split("+")[0] == "aggregate":
			last_snakerule = last_snakerule.split("+")[1]
			for each_input_extension in snakerule_json_Dict[snakerule]["input_extension"]:
				snakerule_input_List.append(OUTPUT + TITLE + "/" + design + "/*/*/output/*" + last_snakerule + "*" + each_input_extension)
			else:
				pass
		else:
			#This is usual rule
			for each_input_extension in snakerule_json_Dict[snakerule]["input_extension"]:
				##
				snakerule_input_List.append(OUTPUT + TITLE + "/" + design + "/*/" + sample + "/output/" + sample + "*" + last_snakerule + "*" + each_input_extension)
			else:
				##
				pass
	return snakerule_input_List


def build_snakerule_target_list(config_Dict, design, sample, snakefile, snakerule, last_snakerule_List):
	"""
	"""
	# ++++++++++++++++++++++++++++++
	TITLE = config_Dict["TITLE"]
	OUTPUT = config_Dict["OUTPUT"]
	# ++++++++++++++++++++++++++++++
	snakerule_json_Dict = get_app_json_dict(config_Dict, snakerule)

	target = OUTPUT + TITLE + "/" + design + "/" + snakefile + "/" + sample + "/target/" + sample + snakerule_json_Dict[snakerule]["target"]
	return target

def build_environment_snakerule_target_list(config_Dict, sample, snakefile, snakerule, last_snakerule_List):
	"""
	"""
	# ++++++++++++++++++++++++++++++
	TITLE = config_Dict["TITLE"]
	OUTPUT = config_Dict["OUTPUT"]
	# ++++++++++++++++++++++++++++++
	snakerule_json_Dict = get_app_json_dict(config_Dict, snakerule)

	target = OUTPUT + TITLE + "/" + snakerule + "/target/" + sample + snakerule_json_Dict[snakerule]["target"]
	return target


def build_snakerule_dict(config_Dict, snakefile, snakerule, last_snakerule_List):
	"""
	"""
	snakerule_Dict = {}
	metadata_file_Path = config_Dict["METADATA"]
	metadata_Dict = build_metadata_dict(metadata_file_Path)
	TITLE = config_Dict["TITLE"]
	End_point_List = []
	# ++++++++++++++++++++++++++++++
	if snakefile in ["build_environment"]:
		##
		snakerule_Dict[snakefile] = {}
		snakerule_Dict[snakefile][TITLE] = {}
		snakerule_Dict[snakefile][TITLE]["input"] = []
		snakerule_Dict[snakefile][TITLE]["target"] = build_environment_snakerule_target_list(config_Dict, TITLE, snakefile, snakerule, last_snakerule_List)
		End_point_List.append(snakerule_Dict[snakefile][TITLE]["target"])
	else:
		##

		for each_design in metadata_Dict:
			##
			snakerule_Dict[each_design] = {}

			if snakefile in ["report"]:
				#
				snakerule_Dict[each_design][TITLE] = {}
				snakerule_Dict[each_design][TITLE]["input"] = build_snakerule_input_list(config_Dict, each_design, TITLE, snakefile, snakerule, last_snakerule_List)
				snakerule_Dict[each_design][TITLE]["target"] = build_snakerule_target_list(config_Dict, each_design, TITLE, snakefile, snakerule, last_snakerule_List)
				End_point_List.append(snakerule_Dict[each_design][TITLE]["target"])
			
			else:
				for each_sample in metadata_Dict[each_design]:
					##
					snakerule_Dict[each_design][each_sample] = {}
					snakerule_Dict[each_design][each_sample]["input"] = build_snakerule_input_list(config_Dict, each_design, each_sample, snakefile, snakerule, last_snakerule_List)
					snakerule_Dict[each_design][each_sample]["target"] = build_snakerule_target_list(config_Dict, each_design, each_sample, snakefile, snakerule, last_snakerule_List)

					End_point_List.append(snakerule_Dict[each_design][each_sample]["target"])
				else:
					##
					pass
		else:
			pass
	return (snakerule_Dict, End_point_List)


def build_snakemake_dict(config_Dict):
	"""
	"""
	snakemake_Dict = {}
	# ++++++++++++++++++++++++++++++++++++
	snakefile_order_List = config_Dict["EXECUTION"]["execution_order"]
	# ++++++++++++++++++++++++++++++++++++
	for snakefile in snakefile_order_List:
		##
		snakerule_order_List = config_Dict["EXECUTION"]["application_order"][snakefile]
		snakemake_Dict[snakefile] = {}
		snakemake_Dict[snakefile]["End_point_List"] = []
		for snakerule in snakerule_order_List:
			##
			#snakerule_json_Dict = get_app_json_dict(config_Dict, snakerule)
			last_snakerule_List = config_Dict["EXECUTION"]["application_dependency"][snakefile][snakerule]
			snakemake_Dict[snakefile][snakerule], End_point_List  = build_snakerule_dict(config_Dict, snakefile, snakerule, last_snakerule_List)
			snakemake_Dict[snakefile]["End_point_List"].extend(End_point_List)
		else:
			##
			pass
	else:
		##
		pass
	return snakemake_Dict

# ################################### SNAKEMAKE EXECUTION ########################### #	

def build_initialize_script(general_Dict):
	"""
	"""
	# ++++++++++++++++++++++++++++++++++++
	initialize_script_String = """#!/bin/bash
		# +++++++++++++++++++++++++++++++++++++++++++
		mkdir -p {MAIN_PATH}
		mkdir -p {REPORT_PATH}
		mkdir -p {LOG_PATH}
		mkdir -p {TEMP_PATH}
		mkdir -p {TARGET_PATH}
		mkdir -p {OUTPUT_PATH}
		mkdir -p {SCRIPT_PATH}
		# -------------------------------------------
	""".format(**general_Dict)
	# ------------------------------------

	return initialize_script_String


def build_finalize_script(argument_Dict, initialize_script_String, main_execution_script):
	"""
	"""

	if main_execution_script is None:
		return ""
	finalize_script_String = ""
	finalize_script_String += initialize_script_String 
	if DEBUG_MODE is True:
		main_execution_script = main_execution_script.replace("\t", "#")
	else:
		pass
	
	finalize_script_String += main_execution_script
	
	return finalize_script_String.replace("\t", "")


def build_general_dict(config_Dict, wildcards):
	"""
	general dict will simplify the requirement of application execution,
	it should be the only input to any application execution 
	"""
	general_Dict = {}
	# ++++++++++++++++++++++++++++++++++++
	#wildcards info
	TITLE = config_Dict["TITLE"]
	SNAKERULE = wildcards.snakerule
	#APP_DICT = config_Dict['APP_DICT']
	SAMPLE = wildcards.sample
	PREFIX = wildcards.prefix
	SNAKEFILE = config_Dict["SNAKEFILE"]
	if SNAKEFILE in ["build_environment"]:
		DESIGN = SNAKEFILE
	else:
		DESIGN = PREFIX.split("/" + SNAKEFILE)[0].split(TITLE + "/")[1]
	# ++++++++++++++++++++++++++++++++++++
	general_Dict["SAMPLE"] = SAMPLE
	general_Dict["TITLE"] = TITLE
	general_Dict["DESIGN"] = DESIGN
	general_Dict["PREFIX"] = PREFIX
	general_Dict["SNAKEFILE"] = SNAKEFILE
	general_Dict["SNAKERULE"] = SNAKERULE
	general_Dict["WORKFLOW_BASEDIR"] = config_Dict["WORKFLOW_BASEDIR"]
	# general_Dict["LAST_SNAKERULE"] = LAST_SNAKERULE
	# ++++++++++++++++++++++++++++++++++++
	#path info
	general_Dict["MAIN_PATH"] = general_Dict["PREFIX"].split("target/")[0]
	general_Dict["REPORT_PATH"] = general_Dict["MAIN_PATH"] + "report/"
	general_Dict["LOG_PATH"] = general_Dict["MAIN_PATH"] + "log/"
	general_Dict["TARGET_PATH"] = general_Dict["MAIN_PATH"] + "target/"
	general_Dict["TEMP_PATH"] = general_Dict["MAIN_PATH"]  + "temp/"
	general_Dict["OUTPUT_PATH"] = general_Dict["MAIN_PATH"] + "output/"
	general_Dict["SCRIPT_PATH"] = general_Dict["MAIN_PATH"] + "script/"
	# ++++++++++++++++++++++++++++++++++++
	#file info
	general_Dict["GENERAL_TAG"] = general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"]
	general_Dict["TARGET_FILE"] = general_Dict["TARGET_PATH"] + general_Dict["SAMPLE"] + "." + general_Dict["SNAKERULE"] + ".target"
	general_Dict["LOG_FILE"] = general_Dict["LOG_PATH"] + general_Dict["GENERAL_TAG"] + ".log"
	general_Dict["SCRIPT_FILE"] = general_Dict["SCRIPT_PATH"] + general_Dict["GENERAL_TAG"] + ".sh"
	
	# ++++++++++++++++
	#metadata
	general_Dict["CONFIG"] = config_Dict
	# ++++++++++++++++
	#I/O
	item_List = []
	snakerule_input_List = config_Dict["IO"][SNAKEFILE][SNAKERULE][DESIGN][SAMPLE]["input"]
	
	for each_input in snakerule_input_List:
		##
		item_List.extend(glob.glob(each_input))
	else:
		pass
	item_List = list(set(item_List))
	item_List.sort()
	general_Dict["INPUT"] = item_List
	# ++++++++++++++++++++++++++
	return general_Dict


def build_execution_dict(config_Dict, wildcards):
	"""
	"""
	general_Dict = build_general_dict(config_Dict, wildcards)

	add_date_tag(general_Dict["LOG_FILE"])
	add_date_tag(general_Dict["SCRIPT_FILE"])

	
	# ++++++++++++++++++++++++++++++++++++++++++
	
	snakerule = wildcards.snakerule
	apps_path = os.path.abspath(general_Dict["WORKFLOW_BASEDIR"] + "/../applications")
	application_name = snakerule.split("_")[0]
	application_path = apps_path + "/" + application_name + "/" + snakerule + ".py"
	application_spec = importlib.util.spec_from_file_location(application_name, application_path)
	application_module = importlib.util.module_from_spec(application_spec)
	application_spec.loader.exec_module(application_module)
	build_execution_function_string = "build_" + snakerule + "_execution_script"
	
	main_execution_script_List = getattr(application_module, build_execution_function_string)(general_Dict)
	
	# +++++++++++++++++++++++++++++++++++++++++++
	
	initialize_script_String = build_initialize_script(general_Dict)

	finalize_script_String = build_finalize_script(general_Dict, initialize_script_String, main_execution_script_List)

	execution_Dict = {}
	execution_Dict["execution_script"] = finalize_script_String
	execution_Dict["general_Dict"] = general_Dict

	return execution_Dict