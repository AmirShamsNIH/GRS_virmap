# ################################### INFO ######################################## #
# Author: Amir Shams
# Project: GRS_VIRmap
# Aim: Snakemake pipeline for viral variation analysis
# Goal: will process/prepare input files for steps downstream
# ################################### IMPORT ###################################### #
library_path = os.path.abspath(workflow.basedir + "/../library")
sys.path.append(library_path)
import snakemake_factory
# ################################### PIPELINE FLOW ############################### #
# ################################### PIPELINE FUNCTION ########################### #

def get_input(wildcards):
	"""
	provide files for snakemake input 
	1- if first snakerule return actual initial input file
	2- if not then provide target file of last_snakerule as input(for priority listing)
	"""
	# +++++++++++++++++++++++++++++++
	title = config["TITLE"]
	snakefile = config["SNAKEFILE"]
	# last_snakefile = config["snakefile"]["last_snakefile"]
	snakerule = wildcards.snakerule
	prefix = wildcards.prefix
	sample = wildcards.sample
	if snakefile in ["build_environment"]:
		design = snakefile
	else:
		design = prefix.split("/" + snakefile)[0].split(title + "/")[1]
	# +++++++++++++++++++++++++++++++
	snakerule_index = app_List.index(snakerule)
	
	if snakerule_index == 0:
		#first
		input_List = []
		for each_input in config["IO"][snakefile][snakerule][design][sample]["input"]:
			##
			input_List.extend(glob.glob(each_input))
		else:
			##
			pass
	else:
		#not first snakerule
		last_snakerule = app_List[snakerule_index - 1]
		input_List =  config["IO"][snakefile][last_snakerule][design][sample]["target"]
	return input_List

# ################################### PIPELINE RULES ############################## #


rule snakemake:
	"""
	"""
	input:
		get_input
	output:
		"{prefix}{sample}.{snakerule}.target"
	wildcard_constraints:
		prefix = "(.+\\/)+",
		# report_path = "(.+\\/)+",
		# log_path = "(.+\\/)+",
		# sample = "((?!.*\\/).*)+",
		# suffix = "((\\.).*)+",
	params:
	resources:
	message: "{wildcards.prefix} | {wildcards.sample} | {wildcards.snakerule} "
	run:
		# ++++++++++++++++++++++++++++++++++++
		execution_Dict = snakemake_factory.build_execution_dict(config, wildcards)
		execution_script = execution_Dict["execution_script"]
		general_Dict = execution_Dict["general_Dict"]
		if execution_script == "":
			print("buidling script was unsuccessful")
		else:
			shell("""
				mkdir -p {general_Dict[SCRIPT_PATH]}
				mkdir -p {general_Dict[TARGET_PATH]}
				mkdir -p {general_Dict[LOG_PATH]}

				printf "%s\\n" '{execution_script}' 2>&1 | tee --append {general_Dict[TARGET_FILE]}.tmp {general_Dict[SCRIPT_FILE]} >/dev/null
				chmod 755 {general_Dict[SCRIPT_FILE]}
			""")
			shell("""
			cd {general_Dict[SCRIPT_PATH]}
			
			printf "%s\\n" 'EXECUTING....' 2>&1 | tee --append {general_Dict[LOG_FILE]} > /dev/null
			
			{general_Dict[SCRIPT_FILE]}
			
			if [ $? -eq 0 ]; then
				mv {general_Dict[TARGET_FILE]}.tmp {general_Dict[TARGET_FILE]}
				printf "%s\\n" 'DONE!!' 2>&1 | tee --append {general_Dict[LOG_FILE]} >/dev/null
			else
				printf "%s\\n" 'FAIL!!' 2>&1 | tee --append {general_Dict[LOG_FILE]} >/dev/null
				echo "Can not move temporary target to permanent target"
			fi
			""")
		# ------------------------------------
# ################################### FINITO ###################################### #