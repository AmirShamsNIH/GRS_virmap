# ################################### INFO ######################################## #
# Author: Amir Shams
# Project: GRS_VIRmap
# Aim: Snakemake pipeline for viral variation analysis
# Goal: will process/prepare input files for steps downstream
# ################################### IMPORT ###################################### #

import os
import sys
import glob
import json
import pathlib

apps_path = os.path.abspath(workflow.basedir + "/../applications")
library_path = os.path.abspath(workflow.basedir + "/../library")
execution_path = os.path.abspath(workflow.basedir + "/../execution")
scripts_path = os.path.abspath(workflow.basedir + "/../scripts")

sys.path.append(apps_path)
sys.path.append(library_path)
sys.path.append(execution_path)
sys.path.append(scripts_path)
import snakemake_factory
# ################################### INCLUDE ##################################### #

# ++++++++++++++++++++++++++++++++++++
execution_file = execution_path + "/snakemake_execution.json"
with open(execution_file) as f:
	execution_Dict = json.load(f)
# ++++++++++++++++++++++++++++++++++++
config["WORKFLOW_BASEDIR"] = workflow.basedir
config["SNAKEFILE"] = "report"
snakefile = "report"
config["EXECUTION"] = execution_Dict
config["IO"] = snakemake_factory.build_snakemake_dict(config)
app_List = execution_Dict["application_order"][snakefile]
# ################################### WILDCARDS ################################### #
# ################################### CONFIGURATION ############################### #
# ################################### CONFIG ###################################### #
# ################################### PIPELINE FLOW ############################### #
rule End_Point:
	input:
		config["IO"][snakefile]["End_point_List"]
# ################################### PIPELINE FUNCTION ########################### #
# ################################### PIPELINE RULES ############################## #
include: "snakemake.smk"
# ################################### FINITO ###################################### #