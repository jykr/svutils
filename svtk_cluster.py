"""
Preprocess MANTA vcf output and run svtk vcfcluster to merge SV calls.


"""

import sys
import os
import subprocess

#if len(sys.argv) != 5: 
#	print("Usage: svtk_cluster.py file_with_manta_vcf_paths merged_record_prefix sample_bnd_ids std_vcf_output_path")
#	exit(1)

# For converting MANTA vcf with inversions as BND into INV
ref_fa = "/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh37d5/human_g1k_v37_decoy.fasta"
manta_inversion_converter = "/home/jr351/projects/Parklab_rotation/SV/Manta/manta-1.6.0.centos6_x86_64/libexec/convertInversion.py"
samtools_path = "/home/jr351/bin/samtools"

"""
vcf_files = open(sys.argv[1], "r").readlines()
merged_record_prefix = sys.argv[2]
sample_bnd_prefix_list = sys.argv[3].strip().split(',')
std_vcf_output_path = sys.argv[4]
restart_file_conversion = True
output_vcf_dir = sys.argv[5]
"""
def convert_inv(vcf_file, output_file, force_do):
	"""Convert BNDs in MANTA output into INVs when possible. Done by script provided in MANTA.
    
    Args:
        vcf_file (str): MANTA output vcf file that contains inversions that are annotated as "BND"s and can be converted into "INV"s. 
        output_file (str): The path of output vcf file where the inversions are converted into "INV"s. 
        force_do (bool): If False, don't re-do conversion where output_file already exists.
    """

	if os.path.isfile(output_file): 
		if not force_do: return
	if "svaba.sv" in vcf_file: 
		print("vcf file is output of svaba, not INV_converting")
		subprocess.run(["cp", vcf_file, output_file])
		return
	conversion_command = manta_inversion_converter + " " + samtools_path + " " + ref_fa + " " + vcf_file + " > " + output_file
	print(conversion_command)
	os.system(conversion_command)
	print("INV converted vcf at: " + output_file)


def sort_vcf(input_vcf, output_vcf):
	subprocess.run(['bcftools', 'sort', input_vcf, '-o', output_vcf + "tmp"])
	subprocess.run(['mv', output_vcf + "tmp", output_vcf])
	return
	


def standardize_vcf(bnd_prefix, vcf_file, output_file, force_do, call_null = False):
	"""Standardize vcf with svtk so it can be clustered by svtk cluster.
    
    Args:
        bnd_prefix (str): Prefix of vcf file IDs in standardized vcf
        vcf_file (str): input vcf file to be svtk standardized
        output_file (str): The path of output vcf file
        force_do (bool): If False, don't re-run svtk standardize if output_file already exists.
        call_null (bool): Include 0/0 variants in the output
    """

	if os.path.isfile(output_file): 
		if not force_do: return
	prefix_arg = ""

	caller = "manta"
	if "svaba.sv" in vcf_file: caller = "svaba"

	if bnd_prefix != "":
		prefix_arg = "-p " + bnd_prefix
		subprocess.run(['svtk', 'standardize', "-p", bnd_prefix, vcf_file, output_file, "manta"])	# svaba output can be parsed by source=manta in svtk standardize
		if caller == "svaba": sort_vcf(output_file, output_file)
		return
	elif call_null: 
		subprocess.run(['svtk', 'standardize',"--call-null-sites", "include-reference-sites", vcf_file, output_file, "manta"])
		if caller == "svaba": sort_vcf(output_file, output_file)	
	else: 
		subprocess.run(['svtk', 'standardize', vcf_file, output_file, "manta"])
		if caller == "svaba": sort_vcf(output_file, output_file)
		
	print("standardized to " + output_file)


def preprocess_vcf_files(output_vcf_dir, vcf_files, sample_bnd_prefix_list, std_vcf_output_path, force_do = True, invert = True):
	""" Preprocess vcf files in vcf_files file.
    0. Convert to "INV" annotation from original MANTA output if invert == True
	1. Standardize with sample_bnd_prefix (.std2.vcf)
	2. Standardize with original bnd IDs (.std.vcf)
	3. Write files in 1 to "mergelist" file

    Args:
        output_vcf_dir (str): directory where the output vcf file will be written
        vcf_files (str): The name of text file containing input vcf file path, one per line.
        sample_bnd_prefix_list (str): The str list of prefixes of vcf IDs in the output standardized vcf files.
        std_vcf_output_path (str): Directory where standardized vcf will be written
        force_do (bool): If False, files won't be generated again if the output file already exsits.
        invert (bool): If False, do not convert inversion annotation from "BND"s to "INV"s (refer to MANTA User Documentation for more detail).
	"""
	std_vcf_list_file = open(output_vcf_dir + "/mergelist", "w")

	for i in range(0, len(vcf_files)):
		vcf_files[i] = vcf_files[i].strip()
		pathsplit = vcf_files[i].split('/')
		sample_id = pathsplit[len(pathsplit)-4]

		# Convert inversions from BND to INV
		vcf_prefix = pathsplit[-1].split('.vcf')[0]
		
		if invert: 	
			vcf_inv_converted = vcf_prefix + ".inv.vcf"
			inversion_converted_vcf_path = "/".join(pathsplit[:-1]) + "/" + vcf_inv_converted
			
			convert_inv(vcf_files[i], inversion_converted_vcf_path, force_do)

			os.system("cp " + inversion_converted_vcf_path + " ./" + vcf_inv_converted)
		else: 
			vcf_inv_converted = pathsplit[-1]
			os.system("cp " + vcf_files[i] + " ./" + vcf_inv_converted)
			if ".gz" in vcf_inv_converted: 
				os.system("gunzip " + vcf_inv_converted)
				vcf_inv_converted = vcf_inv_converted.split('.gz')[0]
		
		call_null = False
		if vcf_prefix == "allSV": call_null = True
		standardize_vcf(sample_bnd_prefix_list[i], vcf_inv_converted, std_vcf_output_path + "/" + sample_id + "." + vcf_prefix + ".std2.vcf", force_do, call_null)
		standardize_vcf("", vcf_inv_converted, std_vcf_output_path + "/" + sample_id + "." +  vcf_prefix + ".std.vcf", force_do)
		
		std_vcf_list_file.write(std_vcf_output_path + "/" + sample_id + "." + vcf_prefix + ".std2.vcf\n")
		subprocess.run(["rm", vcf_inv_converted])

def run_vcfcluster(mergelist_file, merged_record_prefix, svtk_cluster_output_file):
	"""Cluster standardized vcf files"""
	arg_list = ['svtk', 'vcfcluster', '--preserve-ids', "--svtypes", "INS,DEL,DUP,INV,BND", "-d 100", "-f 0.7", "-p", merged_record_prefix, mergelist_file, svtk_cluster_output_file]
	print(arg_list)
	subprocess.run(arg_list)
#command = "svtk vcfcluster --preserve-ids --svtypes 'INS,DEL,DUP,INV,BND' -d 100 -f 0.7 -p "+ merged_record_prefix + " " + std_vcf_output_path + "/mergelist  " + std_vcf_output_path + "/" + merged_record_prefix + ".vcf"

def cluster_manta_output(output_vcf_dir, vcf_list_file, sample_bnd_prefix_list, std_vcf_output_path, clustered_SV_ID_prefix, force_preprocess = True, invert = True):
	"""Main function to call preprocess_vcf_files and run_vcfcluster"""
	vcf_files = open(vcf_list_file, "r").readlines()
	preprocess_vcf_files(output_vcf_dir, vcf_files, sample_bnd_prefix_list, std_vcf_output_path, force_preprocess, invert = invert)
	run_vcfcluster(output_vcf_dir + "/mergelist", clustered_SV_ID_prefix, output_vcf_dir + "/" + clustered_SV_ID_prefix + ".vcf")



