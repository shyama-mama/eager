#!/usr/bin/env python3

# Written by Shyamsundar Ravishankar 
#   Australian Centre for Ancient DNA (ACAD) 
#   University of Adelaide 

"""Script to merge metrics files from picard markduplicates and 
dedup 
"""
import sys
import json
import argparse
import os
import textwrap
import pandas as pd 
import json
import statistics

parser = argparse.ArgumentParser(prog='merge_dedup_metrics.py',
   usage='python %(prog)s [-h] -d [markduplicates|dedup] -i <input directory (default: current directory)> -s <suffix of files> -o <output file name (default: out.<suffix>)>',
   formatter_class=argparse.RawDescriptionHelpFormatter,
   description=textwrap.dedent('''\
   author:
     Shyamsundar Ravishankar (shyamsundar.ravishankar@adelaide.edu.au)

   description:
     %(prog)s Script to merge metrics files from picard markduplicates and dedup 
   '''))
parser.add_argument( "-d", "--dedupper", dest="dedupper", default=None, help="Which dedupper the file is from [markduplicates|dedup]", required=True)
parser.add_argument( "-i", "--input_dir", dest="input_dir", default=os.getcwd(), help="Directory containing sharded metrics files. default: current directory")
parser.add_argument( "-s", "--suffix", dest="suffix", default=None, help="Suffix for files to look for. Usually .metrics, .json, .hist", required=True)
parser.add_argument("-o", "--output_file", dest="output_file", default=None, help="Output file to write the merged metrics to. default: out.<suffix>")
parser.add_argument("-l", "--library", dest="library", default="Unknown library", help="Name of library", required=True)
args = parser.parse_args()

dedupper = args.dedupper
library = args.library
if dedupper != "markduplicates" and dedupper != "dedup":
    print("ERROR: dedupper needs to be 'markduplicates' or 'dedup', got " + dedupper + " instead")
    parser.print_usage()
    sys.exit(1)

suffix = args.suffix 

if not args.input_dir:
    input_dir = os.getcwd()
else:
    input_dir = os.path.abspath(args.input_dir)
    if not os.path.isdir(input_dir):
        print("ERROR: Input directory is not a directory")
        parser.print_usage()
        sys.exit(1)

if not args.output_file:
    output_file = "out" + str(suffix)
else:
    output_file = args.output_file

# Get files with suffix in input_dir 
files_with_suffix = []
# Iterate directory
for file in os.listdir(input_dir):
    # check only text files
    if file.endswith(suffix):
        files_with_suffix.append(file)

no_files = len(files_with_suffix)
print("INFO: Found " + str(no_files) + " files  with suffix " + str(suffix) + " in dir " + input_dir)


#### For Markduplicate metrics read file and combine into a data frame 
# For files in file list get the metrics line into a dataframe 
if dedupper == "markduplicates":
    # Markdup columns
    columns = ["LIBRARY",
        "UNPAIRED_READS_EXAMINED",
        "READ_PAIRS_EXAMINED",
        "SECONDARY_OR_SUPPLEMENTARY_RDS",
        "UNMAPPED_READS",
        "UNPAIRED_READ_DUPLICATES",
        "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES", 
        "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE"]

    # Create df with columns     
    df = pd.DataFrame(columns=columns)

    # Tracker for fist file to get the first few lines
    first_file = True
    for curr_file in files_with_suffix:
        with open(curr_file, 'r') as f, open(output_file, 'a+') as o:
            found = 0
            for line in f:
                if found == 1:
                    new_row = line.strip().split('\t')
                    new_row.append("NA")
                    new_df = pd.DataFrame([new_row], columns=columns)
                    df = pd.concat([df, new_df])
                    break
                else:
                    if first_file:
                        o.write(line)
                if line.startswith("LIBRARY"):                
                    found = 1
        first_file = False

    # Merge the metrics in the df 
    merged_metrics = []
    for col in columns:
        if col == "LIBRARY":
            merged_metrics.append(library)
        elif col == "PERCENT_DUPLICATION":
            # May potentially have to change this - this doesn't seem right 
            if int(merged_metrics[1]) == 0:
                percent_dup = 0
            else:
                percent_dup = int(merged_metrics[5])/int(merged_metrics[1])
            merged_metrics.append(percent_dup)
        elif col == "ESTIMATED_LIBRARY_SIZE":
            merged_metrics.append("")
        else:
            df[col] = pd.to_numeric(df[col])
            merged_metrics.append(df[col].sum())

    # Write merged metrics to output file
    with open(output_file, 'a') as o:
        o.write('\t'.join(str(e) for e in merged_metrics) + '\n')    
elif dedupper == "dedup":
    # Check if hist or json 
    if "hist" in suffix:
        cummulativeHist = {}
        for curr_file in files_with_suffix:
            with open(curr_file, 'r') as f:
                for line in f:
                    list_line = line.split('\t')
                    key = str(list_line[0])
                    value = int(list_line[1])
                    if key in cummulativeHist:
                        cummulativeHist[key] = cummulativeHist[key] + value
                    else:
                        cummulativeHist[key] = value 
        with open(output_file, 'w') as f: 
            for key, value in cummulativeHist.items(): 
                f.write('%s\t%s\n' % (key, value))
    elif "json" in suffix:
        cumulativejson = {}
        cumulativeFactor = []
        first_file = True
        for curr_file in files_with_suffix:
            with open(curr_file, 'r') as f:
                data = json.load(f)
                if first_file:
                    cumulativejson = data
                    cumulativejson['metadata']['sample_name'] = library
                    first_file = False
                else:
                    cumulativejson['metrics']['total_reads'] = int(cumulativejson['metrics']['total_reads']) + int(data['metrics']['total_reads'])
                    cumulativejson['metrics']['mapped_reads'] = int(cumulativejson['metrics']['mapped_reads']) + int(data['metrics']['mapped_reads'])
                    cumulativejson['metrics']['merged_removed'] = int(cumulativejson['metrics']['merged_removed']) + int(data['metrics']['merged_removed'])
                    cumulativejson['metrics']['forward_removed'] = int(cumulativejson['metrics']['forward_removed']) + int(data['metrics']['forward_removed'])
                    cumulativejson['metrics']['total_removed'] = int(cumulativejson['metrics']['total_removed']) + int(data['metrics']['total_removed'])
                    cumulativejson['metrics']['reverse_removed'] = int(cumulativejson['metrics']['reverse_removed']) + int(data['metrics']['reverse_removed'])
                cumulativeFactor.append(float(data['metrics']['clusterfactor']))
        if cumulativejson['metrics']['total_reads'] == 0:
            cumulativejson['metrics']['dup_rate'] = 0
        else:
            cumulativejson['metrics']['dup_rate'] = "{:.8f}".format(float(cumulativejson['metrics']['total_removed'] / cumulativejson['metrics']['total_reads']))
        cumulativejson['metrics']['clusterfactor'] = str(statistics.median(cumulativeFactor))
        json_object = json.dumps(cumulativejson, indent = 4) 
        with open(output_file, 'w') as f: 
            f.write(json_object)
