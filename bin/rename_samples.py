import pandas as pd
import argparse
import glob
import os
import shutil
import re 

def create_rename_dictionary(df, column_name):

	df2 = df.loc[df[column_name].str.startswith("CFIA")].copy()
	df2['old_name'] = df2[column_name].str.extract('-\d(\d{3})-').squeeze().astype(str)
	df2['old_name'] = df2['old_name'].map(lambda x : x[1:] if x.startswith('0') else x )

	df2['old_name'] = 'C' + df2['old_name']

	rename_dict = df2.set_index('old_name')[column_name].to_dict()

	return rename_dict

def main(metadata_path, directory_name, id_column, execute):
	# Read the CSV metadata sheet using Pandas
	df = pd.read_csv(metadata_path)
	
	# Let's assume the column you want to load into the rename dictionary is named 'filename'.
	# You can change 'filename' to whatever your column name is.
	rename_dict = create_rename_dictionary(df, id_column)

	# Use glob.glob to get all files in the given directory
	files = glob.glob(os.path.join(directory_name, "*"))

	for file_path in files:
		base_file_name = os.path.basename(file_path)

		isolate_name, suffix = base_file_name.split("_S")

		# Check if the file is in the rename_dictionary
		if isolate_name in rename_dict:
			new_isolate_name = rename_dict[isolate_name]
			new_file_path = os.path.join(directory_name, new_isolate_name + "_S" + suffix)
			# Rename the file
			if execute:
				shutil.move(file_path, new_file_path)
			print(f"Renamed {file_path} to {new_file_path}")

	if not execute:
		print("Dry run only. No renaming performed.")


def init_parser():
	parser = argparse.ArgumentParser(description="Rename files in a directory based on a metadata sheet.")
	parser.add_argument("-m", "--metadata", help="Path to the CSV metadata sheet", required=True)
	parser.add_argument("-i", "--id-column", help="Name of the ID column in the metadata sheet", default='Name')
	parser.add_argument("-d", "--directory", help="Path to the directory containing the files to rename", required=True)
	parser.add_argument("-x", "--execute", help="Indicates whether to execute the renaming. Otherwise dry run is shown", action='store_true')
	return parser

if __name__ == "__main__":

	parser = init_parser()
	args = parser.parse_args()

	main(args.metadata, args.directory, args.id_column, args.execute)
