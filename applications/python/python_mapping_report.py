import os
import sys
import glob
import pandas

"""
python /data/shamsaddinisha/Development/GRS_virmap/library/collect_report.py ./test_virmap/Group1 ON563414.3-OP470374.1 viral_phix
"""

def parse_samtools_flagstat_report(samtools_flagstat_report):
	"""
	"""
	with open(samtools_flagstat_report) as f:
		first_line = f.readline().strip('\n')
	mapped_read = first_line.split(" +")[0]
	return mapped_read

def build_mapping_dict(mapping_Dict, output_directory, target_genome_string):
	"""
	"""
	target_genome_List= target_genome_string.split("-")
	for each_genome in target_genome_List:
		##
		if each_genome not in mapping_Dict:
			#
			mapping_Dict[each_genome] = {}
		the_path = output_directory + "/**/report/*" + each_genome + ".samtools_flagstat.txt"
		samtools_report_List = glob.glob(the_path,  recursive=True)
		for each_sample_report in samtools_report_List:

			samtools_report_basename = os.path.basename(each_sample_report).split(".")[0]
			mapped_read = parse_samtools_flagstat_report(each_sample_report)
			mapping_Dict[each_genome][samtools_report_basename] = mapped_read
	else:
		pass

	#print(mapping_Dict)
	return mapping_Dict


def build_contamination_dict(mapping_Dict, output_directory, contamination_genome_string):
	"""
	"""
	contamination_genome_List= contamination_genome_string.split("-")
	for each_genome in contamination_genome_List:
		##
		if each_genome not in mapping_Dict:
			#
			mapping_Dict[each_genome] = {}
		the_path = output_directory + "/**/report/*" + each_genome + ".samtools_flagstat.txt"
		samtools_report_List = glob.glob(the_path,  recursive=True)
		for each_sample_report in samtools_report_List:
			samtools_report_basename = os.path.basename(each_sample_report).split(".")[0]
			mapped_read = parse_samtools_flagstat_report(each_sample_report)
			mapping_Dict[each_genome][samtools_report_basename] = mapped_read
	else:
		pass

	print(mapping_Dict)
	return mapping_Dict


def convert_to_csv(output_directory, mapping_Dict):
	"""
	"""
	# for each_genome in mapping_Dict["Target"]:

	dataframe = pandas.DataFrame.from_dict(mapping_Dict, orient='index').T
	# dataframe.index.name = each_genome
	#print(dataframe)
	dataframe.to_csv(output_directory + "/report/aggregate/output/mapping_report.txt", sep="\t")
	return True




mapping_Dict = {}
output_directory = sys.argv[1]
target_genome_string = sys.argv[2]
contamination_genome_string = sys.argv[3]

mapping_Dict = build_mapping_dict(mapping_Dict, output_directory, target_genome_string)
mapping_Dict = build_contamination_dict(mapping_Dict, output_directory, contamination_genome_string)

convert_to_csv(output_directory, mapping_Dict)
print("Done.")