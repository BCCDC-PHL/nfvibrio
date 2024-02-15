#%%
import os 
import sys
import pandas as pd 
import argparse
from glob import glob 


def init_parser():
	parser = argparse.ArgumentParser()
	# parser.add_argument( 'metadata',  help='FASTA consensus sequence containing all query sequences. Segment MUST be included in the header.')	
	parser.add_argument( 'fastq_path',  help='BLAST results from vsearch, containing Genbank seqs as queries against an annotated GISAID DB')	
	parser.add_argument( '-x', '--execute', action='store_true', help='Run the renaming. Without this flag, will just perform dry runs.')
	
	return parser


def rename(path, execute=False):
	os.chdir(path)
	if execute:
		print("Renaming files...")

	folder_names = glob( '*')

	for filepath in glob( "*/*.gz"):
		fields = filepath.split("_")
		folder_name, file_name = fields[0].split("/")
		
		if not file_name.startswith(folder_name):


			newname = "_".join([folder_name, *fields[1:]])
			#print("Starting to rename file")
			print(filepath, " --> ", os.path.join(os.path.dirname(filepath), newname))


			if execute:
				os.rename(filepath, os.path.join(os.path.dirname(filepath), newname))

	if execute:
		print("Renamed all files.")

		print("Complete.")

	else:
		print("DRY RUN ONLY.")
	
def main(args):
	rename(args.fastq_path,  args.execute)
	exit(0)
if __name__ == '__main__':
	parser = init_parser()
	main(parser.parse_args())
	
