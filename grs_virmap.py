#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# ################################### INFO ####################################### #
# Author: Amir Shams
# Project: GRS_VIRMAP
# Aim: pipeline for viral variation analysis
# ################################### IMPORT ##################################### #

import os
import sys
import re
import argparse
import json
import glob
import subprocess
import pandas
# ################################### CONSTANT ################################## #
__author__   = 'Amir Shams'
__version__  = 'v0.0.1'
__email__    = 'shamsaddinisha@nih.gov'
__home__     =  os.path.dirname(os.path.abspath(__file__))
_name       = os.path.basename(sys.argv[0])
_description = 'pipeline to process viral variation'
# ################################### FUNCTIONS ################################ #
# ################################### FINITO ################################### #

def is_file_exist(file_Path):
	"""
	"""
	
	if os.path.isfile(file_Path) and os.path.exists(file_Path) and os.access(file_Path, os.R_OK):
		return True
	else:
		return False


def is_path_empty(the_Path):
	"""
	"""
	if os.path.exists(the_Path) and os.access(the_Path, os.R_OK) and os.listdir(the_Path):
		#
		return False
	else:
		return True


def is_path_readable(the_Path):
	"""
	"""
	if os.path.exists(the_Path) and os.access(the_Path, os.R_OK):
		#
		return True
	else:
		return False


def is_path_writeable(the_Path):
	"""
	"""
	if os.path.exists(the_Path) and os.access(the_Path, os.W_OK):
		#
		return True
	else:
		return False

def verify_input_path(input_directiry_List):
	"""
	"""
	for each_input in input_directiry_List:
		"""
		"""
		if is_path_readable(each_input) is False:
			print("input path is not accessible; either not exist or not readable!!")
			print("input_path:", each_input)
			print("Aborting!!")
			sys.exit(2)
		else:
			input_List = []
			for each_extension in ["fq", "fastq", "fa", "fasta"]:
				#
				input_List.extend(glob.glob(each_input + "/*." + each_extension + "*"))
			else:
				if len(input_List) == 0:
					print("No fastq/fasta file detected within input_path; please make sure the input_path is correct!")
					print("we support fq, fastq, fa, fasta only")
					print("input_path:", each_input)
					print("Aborting!!")
					sys.exit(2)
				else:
					pass
	return True


def verify_output_path(output_directory):
	"""
	"""
	if is_path_readable(output_directory):
		if is_path_writeable(output_directory):
			pass
		else:
			print("output path is not writable; please check permission")
			print("output_path:", output_directory)
			print("Aborting!!")
			sys.exit(2)
	else:
		print("output path is not accessible; either not exist or not readable!!")
		print("output_path:", output_directory)
		print("Aborting!!")
		sys.exit(2)
	
	return True


def verify_metadata(metadata_file_path):
	"""
	"""
	if is_file_exist(metadata_file_path) is False:
		print("metadata file path is not avaiable or accessible")
		print("metadata:", metadata_file_path)
		print("Aborting!!")
		sys.exit(2)
	
	return True


def verify_execution_mode(execution_mode):
	"""
	"""
	if execution_mode.upper() not in ["SLURM", "LOCAL"]:
		print("execution mode should be SLURM or LOCAL")
		print("MODE:", execution_mode)
		print("Aborting!!")
		sys.exit(2)
	return True

def verify_execution_platform(execution_platfrom):
	"""
	"""
	if execution_platfrom.upper() not in ["BIOWULF", "BIGSKY"]:
		print("execution mode should be BIOWULF or BIGSKY")
		print("PLATFORM:", execution_platfrom)
		print("Aborting!!")
		sys.exit(2)
	return True


def verify_target_reference(target_reference_Dict, target_reference_index_List):
	"""
	"""
	for each_target in target_reference_index_List:
		if each_target not in target_reference_Dict:
			#
			print("target reference is missing")
			print("TARGET_REFERENCE_INDEX:", each_target)
			print("target_reference_file:", target_reference_Dict.keys())
			print("Aborting!!")
			sys.exit(2)
		target_fasta_file_path = target_reference_Dict[each_target]["fasta"]
		
		if is_file_exist(target_fasta_file_path) is False:
			#
			print("target reference file is missing or not accessible")
			print("TARGET_REFERENCE:", each_target)
			print("target_reference_file:", target_fasta_file_path)
			print("Aborting!!")
			sys.exit(2)
		target_gtf_file_path = target_reference_Dict[each_target]["gtf"]
		
		if is_file_exist(target_gtf_file_path) is False:
			#
			print("target reference file is missing or not accessible")
			print("TARGET_REFERENCE:", each_target)
			print("target_reference_file:", target_gtf_file_path)
			print("Aborting!!")
			sys.exit(2)
	else:
		pass
	return True


def verify_target_reference_path(target_reference_path):
	"""
	"""
	# ++++++++
	if is_path_readable(target_reference_path) is False:
		print("target_reference_path is not accessible; either not exist or not readable!!")
		print("target_reference_path:", target_reference_path)
		print("Aborting!!")
		sys.exit(2)
	# ++++++++
	snpEFF_config_file_path = target_reference_path + "/snpEff.config"
	if is_file_exist(snpEFF_config_file_path) is False:
		#
		print("snpEff.config is missing in targete reference path")
		print("snpEFF_config_file_path:", snpEFF_config_file_path)
		print("Aborting!!")
		sys.exit(2)
	return True
	

def verify_target_list(target_reference_path, target_List):
	"""
	"""
	snpEFF_config_file_path = target_reference_path + "/snpEff.config"
	snpEFF_DF = pandas.read_csv(snpEFF_config_file_path, sep=": ", header=None, engine='python')
	for each_target in target_List:
		##
		if each_target not in snpEFF_DF[1].values:
			#
			print("target is not listed in target reference")
			print("target_reference_path:", target_reference_path)
			print("snpEFF_config_file_path:", snpEFF_config_file_path)
			print("target:", each_target)
			print("Aborting!!")
			sys.exit(2)
	else:
		pass
	return True


def verify_contamination_reference(target_reference_Dict, contamination_reference_index_List):
	"""
	"""
	for each_target in contamination_reference_index_List:
		if each_target not in target_reference_Dict:
			#
			print("contamination reference is missing")
			print("CONTAMINATION_REFERENCE_INDEX:", each_target)
			print("contamination_reference:", target_reference_Dict.keys())
			print("Aborting!!")
			sys.exit(2)
		target_fasta_file_path = target_reference_Dict[each_target]
		target_fasta_path = os.path.dirname(target_fasta_file_path)
		if is_path_empty(target_fasta_path):
			#
			print("contamination reference file is missing or not accessible")
			print("CONTAMINATION_REFERENCE:", each_target)
			print("contamination_reference:", target_fasta_file_path)
			print("Aborting!!")
			sys.exit(2)
	else:
		pass
	return True


def verify_kraken2_reference(target_reference_Dict, kraken2_reference_index_List):
	"""
	"""
	for each_target in kraken2_reference_index_List:
		if each_target not in target_reference_Dict:
			#
			print("kraken2 reference is missing")
			print("kraken2_REFERENCE_INDEX:", each_target)
			print("kraken2_reference:", target_reference_Dict.keys())
			print("Aborting!!")
			sys.exit(2)
		target_fasta_file_path = target_reference_Dict[each_target]
		if is_path_empty(target_fasta_file_path):
			#
			print("kraken2 reference file is missing or not accessible")
			print("kraken2_REFERENCE:", each_target)
			print("kraken2_reference_file:", target_fasta_file_path)
			print("Aborting!!")
			sys.exit(2)
	else:
		pass
	return True


def verify_parameters(parameters_Dict):
	"""
	"""
	return True


def validate_input_json(json_file_path):
	"""
	"""
	# ++++++++++++++++++++++++++++
	f = open(json_file_path, "r")
	input_json_Dict = json.load(f)
	# ++++++++++++++++++++++++++++
	if "TITLE" in input_json_Dict:
		pipeline_title = input_json_Dict["TITLE"]
	else:
		print("TITLE value is required!")
		print("Aborting!!")
		sys.exit(2)
	# ++++++++++++++++++++++++++++
	if "INPUT" in input_json_Dict:
		input_directiry_List = input_json_Dict["INPUT"]
		verify_input_path(input_directiry_List)
	else:
		print("INPUT value is required!")
		print("Aborting!!")
		sys.exit(2)
	# ++++++++++++++++++++++++++++
	if "OUTPUT" in input_json_Dict:
		output_directory = input_json_Dict["OUTPUT"]
		verify_output_path(output_directory)
	else:
		print("OUTPUT value is required!")
		print("Aborting!!")
		sys.exit(2)
	# ++++++++++++++++++++++++++++
	if "EXECUTION_MODE" in input_json_Dict:
		execution_mode = input_json_Dict["EXECUTION_MODE"]
		verify_execution_mode(execution_mode)
	else:
		print("EXECUTION_MODE value is required!")
		print("Aborting!!")
		sys.exit(2)
	# ++++++++++++++++++++++++++++
	if "EXECUTION_PLATFORM" in input_json_Dict:
		execution_platform = input_json_Dict["EXECUTION_PLATFORM"]
		verify_execution_platform(execution_platform)
	else:
		print("EXECUTION_PLATFORM value is required!")
		print("Aborting!!")
		sys.exit(2)
	
	# ++++++++++++++++++++++++++++
	if "METADATA" in input_json_Dict:
		metadata_file_path = input_json_Dict["METADATA"]
		verify_metadata(metadata_file_path)
	# ++++++++++++++++++++++++++++
	if "TARGET_LIST" in input_json_Dict:
		#
		target_List = input_json_Dict["TARGET_LIST"]
		if "TARGET_REFERENCE_PATH" in input_json_Dict:
			target_reference_path = input_json_Dict["TARGET_REFERENCE_PATH"]
			verify_target_reference_path(target_reference_path)
			verify_target_list(target_reference_path, target_List)
			#target_reference_Dict = input_json_Dict["TARGET_REFERENCE"][execution_platform]
			#verify_target_reference(target_reference_Dict, target_reference_index_List)
		else:
			print("TARGET_REFERENCE_PATH value is required!")
			print("Aborting!!")
			sys.exit(2)
	else:
		print("TARGET_LIST value is required!")
		print("Aborting!!")
		sys.exit(2)

	# ++++++++++++++++++++++++++++
	print("input_json validated!")

	return True
	

def main(args):
	"""
	"""
	json_file_path = args['input_json']
	
	validate_input_json(json_file_path)

	pipeline_directory = os.path.dirname(os.path.realpath(__file__))
	build_execution_script_python_file_path = pipeline_directory + "/execution/build_execution_script.py"
	process = subprocess.Popen(
		['python', build_execution_script_python_file_path, json_file_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE
		)
	stdout, stderr = process.communicate()
	print(stdout)
	print(stderr)
	if int(process.returncode) == 0:
		print("execution script created succcessfully!")
		print("Please submit the execution script into server cluster!")
	else:
		print("execution script building failed!")

	
	return True

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Usage: grs_virmap.py <input_json>")
	parser.add_argument('-i','--input_json', help='metadata csv file; required for analysis', required=True)
	
	args = vars(parser.parse_args())
	
	if len(args) < 1:
		print("you need to prvide I/O for running pipeline")
		sys.exit(2)
	

	main(args)