# ################################### INFO ######################################## #
# Author: Amir Shams
# Project: a project
# Aim: ...
# ################################### IMPORT ###################################### #


import os
import sys
import json
import glob
import re
# ################################### INCLUDE ##################################### #
# ################################### FUNCTIONS ################################### #
# ################################### CONFIGURATION ############################### #
# ################################### WILDCARDS ################################### #
# ################################### PIPELINE FLOW ############################### #


def build_metadata(input_json_Dict):
	"""
	"""
	fastq_extension_regex = "(\.fastq\.gz$|\.fq\.gz$|\.fastq$|\.fq$)"
	fasta_extension_regex = "(\.fasta\.gz$|\.fa\.gz$|\.fasta$|\.fa$)"
	forward_name_tag_regex = "(_R1|\.r1|\.R1|\.1\.).*"
	reverse_name_tag_regex = "(_R2|\.r2|\.R2|\.2\.).*"

	input_directory_List = input_json_Dict["INPUT"]
	output_directory = input_json_Dict["OUTPUT"]
	TITLE = input_json_Dict["TITLE"]
	metadata_file_path = output_directory + "/" + TITLE + ".sample_metadata.csv"
	input_List = []
	for each_input in input_directory_List:
		"""
		"""
		for each_extension in ["fq", "fastq", "fa", "fasta"]:
			#
			input_List.extend(glob.glob(each_input + "/*." + each_extension + "*"))
		else:
			pass
	else:
		pass
	input_List.sort()
	metadata_Dict = {}
	sample_List = []
	for each_item in input_List:
		##
		item_name = os.path.basename(each_item)
		item_name_no_extension = re.sub(fastq_extension_regex, "", item_name)
		item_name_no_extension = re.sub(fasta_extension_regex, "", item_name_no_extension)
		item_name_no_extension_no_tag = re.sub(forward_name_tag_regex, "", item_name_no_extension)
		item_name_no_extension_no_tag = re.sub(reverse_name_tag_regex, "", item_name_no_extension_no_tag)
		sample_List.append(item_name_no_extension_no_tag)
	else:
		pass
	sample_List = list(set(sample_List))
	metadata_String = "Sample,Design,Type\n"
	for each_sample in sample_List:
		##
		metadata_String += each_sample + ",Group1,CASE\n"
	else:
		pass
	o = open(metadata_file_path, "w")
	o.write(metadata_String)
	o.close()
	return metadata_file_path


def build_execution_snakemake_script(json_file_path, input_json_Dict):
	"""
	"""
	execution_bash_script = ""
	execution_bash_script += """#!/bin/bash
	# +++++++++++++++++++++++++++++++++++
	# sbatch --cpus-per-task=16 --mem=16g --time=24:00:00 {TITLE}.execution_script.sh
	if [[ "{PLATFORM}" == "BIOWULF" ]]
	then
	module load snakemake || exit 1
	snakemake=$(which snakemake) || exit 1
	else
	snakemake=/gs1/RTS/NextGen/bin/grs_virmap/bin/snakemake || exit 1
	fi
	# -----------------------------------
	# +++++++++++++++++++++++++++++++++++
	echo 'GRS_virmap Pipeline execution initiated at: '$(date)
	mkdir -p {WORKDIR} || exit 1
	mkdir -p {WORKDIR}/snakemake_log
	mkdir -p {WORKDIR}/snakemake_log/pre_process
	mkdir -p {WORKDIR}/snakemake_log/alignment
	mkdir -p {WORKDIR}/snakemake_log/variant_caller
	mkdir -p {WORKDIR}/snakemake_log/report
	# -----------------------------------
	# +++++++++++++++++++++++++++++++++++
	cd {WORKDIR}/snakemake_log
	# -----------------------------------
	# +++++++++++++++++++++++++++++++++++
	""".format(
		WORKDIR=input_json_Dict["OUTPUT"] + "/" + input_json_Dict["TITLE"],
		TITLE=input_json_Dict["TITLE"],
		PLATFORM=input_json_Dict["EXECUTION_PLATFORM"],
	)
	# -----------------------------------
	# +++++++++++++++++++++++++++++++++++
	pipeline_directory = os.path.dirname(__file__) + "/.." 
	execution_bash_script += """
	$snakemake \\
	--snakefile {PIPELINE_DIR}/snakemake/pre_process.smk \\
	--configfile {execution_json} \\
	--cores --unlock
	""".format(
		PIPELINE_DIR=pipeline_directory,
		execution_json=json_file_path
	)


	if input_json_Dict["EXECUTION_PLATFORM"] == "BIOWULF":
		#

		cluster_slurm_script = """--cluster="sbatch --cpus-per-task={core} --mem={memory} \\
			--partition={partition} --time={time} --mail-type=FAIL --job-name={jobname} \\
			--output={WORKDIR}/snakemake_log/{output} \\
			--error={WORKDIR}/snakemake_log/{error} {extra}" """.format(
			WORKDIR=input_json_Dict["OUTPUT"] + "/" + input_json_Dict["TITLE"],
			core="{cluster.core}",
			memory="{cluster.memory}",
			partition = "{cluster.partition}",
			time="{cluster.time}",
			jobname="{cluster.jobname}",
			output="{cluster.output}",
			error="{cluster.error}",
			extra = "{cluster.extra}"
		)
	else:
		cluster_slurm_script = """--cluster="sbatch --cpus-per-task={core} --mem={memory} \\
			--time={time} --mail-type=FAIL --job-name={jobname} \\
			--output={WORKDIR}/snakemake_log/{output} \\
			--error={WORKDIR}/snakemake_log/{error}" """.format(
			WORKDIR=input_json_Dict["OUTPUT"] + "/" + input_json_Dict["TITLE"],
			core="{cluster.core}",
			memory="{cluster.memory}",
			time="{cluster.time}",
			jobname="{cluster.jobname}",
			output="{cluster.output}",
			error="{cluster.error}"
		)

	if input_json_Dict["EXECUTION_MODE"] == "local":
		#

		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/build_environment.smk \\
		--configfile {execution_json} \\
		--cores all  --keep-incomplete --rerun-triggers mtime
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		)


		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/pre_process.smk \\
		--configfile {execution_json} \\
		--cores all  --keep-incomplete --rerun-triggers mtime
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		)

		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/alignment.smk \\
		--configfile {execution_json} \\
		--cores all  --keep-incomplete --rerun-triggers mtime
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		)

		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/variant_calling.smk \\
		--configfile {execution_json} \\
		--cores all  --keep-incomplete --rerun-triggers mtime
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		)

		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/report.smk \\
		--configfile {execution_json} \\
		--cores all --keep-incomplete --rerun-triggers mtime
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		)
	
	elif input_json_Dict["EXECUTION_MODE"] == "slurm":


		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/build_environment.smk \\
		--configfile {execution_json} \\
		--cores all  --keep-incomplete --rerun-triggers mtime
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		)
		#
		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/pre_process.smk \\
		--configfile {execution_json} \\
		--cluster-config {PIPELINE_DIR}/execution/cluster_production.yaml \\
		--jobs=10 --max-jobs-per-second=1 --max-status-checks-per-second=0.01 --latency-wait=120 \\
		--keep-incomplete --rerun-triggers mtime  \\
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		) + cluster_slurm_script

		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/alignment.smk \\
		--configfile {execution_json} \\
		--cluster-config {PIPELINE_DIR}/execution/cluster_production.yaml \\
		--jobs=10 --max-jobs-per-second=1 --max-status-checks-per-second=0.01 --latency-wait=120 \\
		--keep-incomplete --rerun-triggers mtime  \\
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		) + cluster_slurm_script

		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/variant_calling.smk \\
		--configfile {execution_json} \\
		--cluster-config {PIPELINE_DIR}/execution/cluster_production.yaml \\
		--jobs=10 --max-jobs-per-second=1 --max-status-checks-per-second=0.01 --latency-wait=120 \\
		--keep-incomplete --rerun-triggers mtime  \\
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		) + cluster_slurm_script

		execution_bash_script += """
		# -----------------------------------
		# +++++++++++++++++++++++++++++++++++
		$snakemake \\
		--snakefile {PIPELINE_DIR}/snakemake/report.smk \\
		--configfile {execution_json} \\
		--cluster-config {PIPELINE_DIR}/execution/cluster_production.yaml \\
		--jobs=10 --max-jobs-per-second=1 --max-status-checks-per-second=0.01 --latency-wait=120 \\
		--keep-incomplete --rerun-triggers mtime  \\
		""".format(
			PIPELINE_DIR=pipeline_directory,
			execution_json=json_file_path
		) + cluster_slurm_script
	else:
		#
		pass
	execution_bash_script += """\n# -----------------------------------
	echo 'GRS_virmap Pipeline execution initiated at: '$(date)"""
	execution_bash_script = execution_bash_script.replace("\t", "")
	return execution_bash_script



	
	
# ################################### MAIN ######################################## #

def main():
	# Sanity check for usage
	
	if len(sys.argv) == 1:
		# Nothing was provided
		print("input json file are missing")
		print("Aborting!!!")
		sys.exit(2)

	input_json_file_Path = sys.argv[1]
	
	# ++++++++++++++++++++++++++++
	f = open(input_json_file_Path, "r")
	input_json_Dict = json.load(f)
	# ++++++++++++++++++++++++++++
	output_directory = input_json_Dict["OUTPUT"]
	TITLE = input_json_Dict["TITLE"]
	json_file_path = output_directory + TITLE + ".execution_metadata.json"
	if "METADATA" not in input_json_Dict:
		metadata_file_path = build_metadata(input_json_Dict)
		input_json_Dict["METADATA"] = metadata_file_path
	else:
		pass
	# ++++++++++++++++++++++++++++
	pipeline_directory = os.path.dirname(os.path.realpath(__file__))
	parameters_json_file_path = pipeline_directory + "/../template/grs_virmap_parameters.json"
	
	f = open(parameters_json_file_path, "r")
	parameters_json_Dict = json.load(f)

	input_json_Dict["PARAMETERS"] = parameters_json_Dict["PARAMETERS"]
	# +++++++++++++++++++++++++++++
	execution_Json = json.dumps(input_json_Dict, indent=4)
	o=open(json_file_path, "w")
	o.write(execution_Json)
	o.close()
	
	
	

	# +++++++++++++++++++++++++++++
	execution_bash_script = build_execution_snakemake_script(json_file_path, input_json_Dict)
	excution_file_path = output_directory + "/" + TITLE + ".execution_script.sh"
	o=open(excution_file_path, "w")
	o.write(execution_bash_script)
	o.close()
	# +++++++++++++++++++++++++++++
	


if __name__ == '__main__':
	main()
# ################################### FINITO ###################################### #