#######################################################################
#######################################################################
# Annovar2VCF
# Developed by Kijin Kim
# skkujin@gmail.com
#
# repository: github.com/KijinKims/Annovar2VCF
#
#######################################################################
#######################################################################

import argparse
import os, sys
import numpy  as np
import pandas as pd
import subprocess

def set_variant_type(row):
	if row["Ref"] == "-":
		if row["Alt"] == "-":
			pass #Error
		else:
			return "Insertion"
	elif row["Alt"] == "-":
		return "Deletion"
	else:
		if len(row["Ref"]) == 1 and len(row["Alt"]) == 1:
			return "SNP"
		else:
			pass #Error

def chr_write(chromosome, file_obj):
	if chromosome.startswith("chr"):
		file_obj.write(chromosome)
	else:
		file_obj.write("chr"+chromosome)

def stack_position(row, file_obj):
	if row["_variant_type"] == "Insertion":
		chr_write(row["Chr"], file_obj)
		file_obj.write('\t'+str(row["Start"])+'\t'+str(row["Start"])+'\n')
	elif row["_variant_type"] == "Deletion":
		chr_write(row["Chr"], file_obj)
		file_obj.write('\t'+str(row["Start"]-1)+'\t'+str(row["Start"]-1)+'\n')

def main():
	#default tab-delimited file
	#input file parsing
	raw_input_file = open(sys.argv[1])
	
	df = pd.read_table(raw_input_file, header=0, dtype={'Chr':object})
	df = df.fillna('.')
	
	#set variant type
	df = df.assign(_variant_type=df.apply(set_variant_type, axis=1))
	
	#make retrieve_seq_from_fasta.pl input file
	ref_fasta_file_location = sys.argv[2]
	rsff_input_file = open("temp_for_rsff.txt","w")
	df.apply(stack_position, axis=1, file_obj=rsff_input_file)
	
	#run retrieve_seq_from_fasta.pl
	run_args=["retrieve_seq_from_fasta.pl", "-format", "tab", "-seqfile", ref_fasta_file_location, "./temp_for_rsff.txt", "-outfile", "./result_rsff.txt", "--verbose"]
	#subprocess.Popen(run_args)
	
	#replace sequence using result
	with open("result_rsff.txt") as result_file:
		result_lines = result_file.read().splitlines()
		result_seq = result_lines[1::2]
	
	result_index = 0
	
	for index, row in df.iterrows():
		if row["_variant_type"] == "Insertion":
			df.loc[index, "Ref"] = result_seq[result_index]
			df.loc[index, "Alt"] = result_seq[result_index] + row["Alt"]
			result_index += 1
		elif row["_variant_type"] == "Deletion":
			df.loc[index, "Start"] = row["Start"] - 1
			df.loc[index, "Ref"] = result_seq[result_index] + row["Ref"]
			df.loc[index, "Alt"] = result_seq[result_index]
			result_index += 1
		else:
			pass
	
		if result_index == len(result_seq):
			pass #error
	
	##export as vcf form
	export_vcf = open(sys.argv[3],"w")
	
	#header
	export_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	column_names = df.columns.values
	ban_columns = ["_variant_type", "Chr", "Start", "End", "Ref", "Alt"]
	additional_column_names = [x for x in column_names if x not in ban_columns]
	export_vcf.write('\t'+'\t'.join(additional_column_names))
	export_vcf.write('\n')
	
	for index, row in df.iterrows():
		CHROM=row["Chr"]
		POS=row["Start"]
		ID="."
		REF=row["Ref"]
		ALT=row["Alt"]
		QUAL=FILTER=INFO="."
		write_line=[CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO]+df.loc[index,additional_column_names].values.tolist()
		write_line=[str(x) for x in write_line]
		export_vcf.write('\t'.join(write_line)+'\n')

if __name__ == "__main__":
	main()
