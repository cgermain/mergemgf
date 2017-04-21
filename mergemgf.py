from pyteomics import mgf
import os
import sys
import argparse

parser = argparse.ArgumentParser(description="Merge 2 MGF files based on X coordinate criteria")
parser.add_argument('ms2_file', action="store", help="Path to the MS2 mgf file.")
parser.add_argument('ms3_file', action="store", help="Path to the MS3 mgf file.")

MZ_CUTOFF = 140

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
	
	merge_result = merge_mgf_files(ms2_file, ms3_file)

	#add _MS3 to the filename
	merged_mgf_filename = ms2_file.split(".")[0].replace("MS2", "MS2_MS3")+".mgf"

	print "\nWriting merged MGF: " + merged_mgf_filename

	mgf.write(merge_result["merged_mgf"], output=merged_mgf_filename)
	print_merge_stats(merge_result)

def merge_mgf_files(ms2_file, ms3_file):
	ms2_count = 0
	ms3_count = 0
	current_count = 0
	merged_count = 0

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
			ms2_count += 1

	print "Reading MS3 file: " + ms3_file
	with mgf.read(ms3_file) as ms3_reader:
		for ms3_temp in ms3_reader:
			ms3_spectrum_list.append(ms3_temp)
			ms3_count+=1

	#Loop through all MS2/MS3 spectra looking for fuzzy matches.
	for ms2_spectrum in ms2_spectrum_list:
		for ms3_index, ms3_spectrum in enumerate(ms3_spectrum_list):
			if compare_spectrums_with_fuzzy_rt(ms2_spectrum, ms3_spectrum):
				merged_xy = merge_xy_arrays(ms2_spectrum, ms3_spectrum)
				ms2_spectrum['m/z array'] = merged_xy[0]
				ms2_spectrum['intensity array'] = merged_xy[1]					
				merged_count += 1
				#remove the element we just found from the list to avoid dupes and save time
				del ms3_spectrum_list[ms3_index]
				break
		merged_mgf.append(ms2_spectrum) #add no matter if it was merged or not
		current_count += 1
		write_progress_bar(current_count, ms2_count)

	return {
		"merged_mgf":merged_mgf,
		"ms2_count":ms2_count,
		"ms3_count":ms3_count,
		"merged_count":merged_count
		}

def compare_spectrums_with_fuzzy_rt(ms2_spectrum, ms3_spectrum):
	return (ms2_spectrum['params']['pepmass'] == ms3_spectrum['params']['pepmass'] and
		ms2_spectrum['params']['charge'] == ms3_spectrum['params']['charge'] and
		abs(float(ms2_spectrum['params']['rtinseconds'])-float(ms3_spectrum['params']['rtinseconds'])) < 5)

def merge_xy_arrays(ms2_spectrum, ms3_spectrum):
	merge_mz = []
	merge_intensity = []
	
	#Take all of the X,Y pairs that have X values lower than MZ_CUTOFF and add them to the merged array
	#Start with ms3 so they remain in sorted order
	for index, x_value in enumerate(ms3_spectrum['m/z array']):
		if (float(x_value) < MZ_CUTOFF):
			merge_mz.append(x_value)
			merge_intensity.append(ms3_spectrum['intensity array'][index])
		else: # skip the rest if over MZ_CUTOFF
			break # can do this since it is ordered

	for index, x_value in enumerate(ms2_spectrum['m/z array']):
		if (float(x_value) >= MZ_CUTOFF):
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

def print_merge_stats(merge_result):
	print "MS2 Count  : " + str(merge_result["ms2_count"])
	print "MS3 Count  : " + str(merge_result["ms3_count"])
	print "Merged Count: " + str(merge_result["merged_count"])

if __name__ == "__main__":
	main()