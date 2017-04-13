from pyteomics import mgf
import os
import sys
import argparse

parser = argparse.ArgumentParser(description="Merge 2 MGF files based on X coordinate criteria")
parser.add_argument('ms2_file', action="store", help="Path to the MS2 mgf file.")
parser.add_argument('ms3_file', action="store", help="Path to the MS3 mgf file.")

def main():
	#TODO run on a folder
	#TODO pull in each mgf file & check for ms2/ms3
	arguments = parser.parse_args()
	if not os.path.isfile(arguments.ms2_file):
		print "Please check that the path to the ms2 file is correct."
		return
	if not os.path.isfile(arguments.ms3_file):
		print "Please check that the path to the ms2 file is correct."
		return
	ms2_file = arguments.ms2_file
	ms3_file = arguments.ms3_file
	merge_spectrums(ms2_file, ms3_file)

def merge_spectrums(ms2_file, ms3_file):
	current_count = 0
	merged_count = 0
	total_count = 0

	#preloading the files into memory
	#ms2 - so we have a total spectra count for the progress bar
	#ms3 - so we don't have to read in repeatedly per ms2 spectra

	merged_mgf = []
	ms3_spectrum_list = []
	ms2_spectrum_list = []

	print "Reading MS2 file: " + ms2_file
	with mgf.read(ms2_file) as ms2_reader:
		for ms2_temp in ms2_reader:
			ms2_spectrum_list.append(ms2_temp)
			total_count += 1

	print "Reading MS3 file: " + ms3_file
	with mgf.read(ms3_file) as ms3_reader:
		for ms3_temp in ms3_reader:
			ms3_spectrum_list.append(ms3_temp)

	for ms2_spectrum in ms2_spectrum_list:
		for ms3_spectrum in ms3_spectrum_list:
			if compare_spectrums(ms2_spectrum, ms3_spectrum):
				merged_xy = merge_xy_arrays(ms2_spectrum, ms3_spectrum)
				#replace the arrays in the ms2 dictionary
				ms2_spectrum['m/z array'] = merged_xy[0]
				ms2_spectrum['intensity array'] = merged_xy[1]
				merged_count += 1
				#no need to continue searching after a match was found
				break
		merged_mgf.append(ms2_spectrum) #add no matter if it was merged or not
		current_count += 1
		write_progress_bar(current_count, total_count)

	#TODO generate output filename & subdirectory
	merged_mgf_filename = "combined_mgfs.mgf"

	print "\nWriting merged MGF: " + merged_mgf_filename
	mgf.write(merged_mgf, output=merged_mgf_filename)
	print "Merged Count: " + str(merged_count)
	print "Total Count: " + str(total_count)

def compare_spectrums(ms2_spectrum, ms3_spectrum):
	if (ms2_spectrum['params']['pepmass'] == ms3_spectrum['params']['pepmass'] and
		ms2_spectrum['params']['charge'] == ms3_spectrum['params']['charge'] and
		ms2_spectrum['params']['rtinseconds'] == ms3_spectrum['params']['rtinseconds']):
		return True
	else:
		return False

def merge_xy_arrays(ms2_spectrum, ms3_spectrum):
	merge_mz = []
	merge_intensity = []
	
	#Take all of the X,Y pairs that have X values lower than 140 and add them to the merged array
	#Start with ms3 so they remain in sorted order
	for index, x_value in enumerate(ms3_spectrum['m/z array']):
		if (float(x_value) < 140):
			merge_mz.append(x_value)
			merge_intensity.append(ms3_spectrum['intensity array'][index])
		else: # skip the rest if over 140
			break # can do this since it is ordered

	for index, x_value in enumerate(ms2_spectrum['m/z array']):
		if (float(x_value) >= 140):
			merge_mz.append(x_value)
			merge_intensity.append(ms2_spectrum['intensity array'][index])

	return [merge_mz, merge_intensity] 

def write_progress_bar(current_count, total_count):
	width = 60
	complete_length = int(round(width * current_count / float(total_count)))

	percent_complete = round(100.0 * current_count / float(total_count), 1)
	bar = '=' * complete_length + '-' * (width - complete_length)

	sys.stdout.write('Merging: [%s] %s%s \r' % (bar, percent_complete, '%'))
	sys.stdout.flush()

if __name__ == "__main__":
	main()