#!/usr/bin/env python3

import re
import argparse
from deduper_ftns import *


def get_args():
    parser = argparse.ArgumentParser(description="dedupper")
    parser.add_argument("-f", "--file", help="Path to SAM file", required=True, type = str)
    parser.add_argument("-p", "--paired", help="State whether SAM file is created from PAIRED-END library (if paired type: paired/if unpaired type: unpaired)", required=False, type = str)
    parser.add_argument("-u", "--umi", help="Path to file with UMIs", required=True, type = str)
    parser.add_argument("-s", "--sorted", help="State whether SAM file is sorted (if sorted type: sorted/if unsorted type: unsorted)", required=True, type = str)
    return parser.parse_args()

#FILEHANDLING
args = get_args()
INPUT_FILE_NAME = args.file
INPUT_FILE = open(INPUT_FILE_NAME, "r")
INPUT_FILE_NAME = INPUT_FILE_NAME.split("/")[-1]            #CREATE OUTPUT FILEHANDLE
INPUT_FILE_NAME = INPUT_FILE_NAME.split(".")[0]
OUTPUT_FILE = open(INPUT_FILE_NAME + "_dedupped.sam","w")
LIB_UMI_LIST_FILE = open(args.umi, "r")
PAIRED = args.paired
SORTED = args.sorted

#ERROR MESSAGE FOR PAIRED END DATA AND UNSORTED DATA
if PAIRED != None:
    PAIRED = PAIRED.lower()
if PAIRED != None:
    if PAIRED == "paired":
        print("***ERROR***" + "\n" + "Sorry, this program does NOT work for paired-end libraries.")
        exit()

if SORTED != None:
    SORTED = SORTED.lower()
if SORTED != None:
    if SORTED == "unsorted":
        print("***ERROR***" + "\n" + "Sorry, this program only works with sorted SAM files.")
        exit()

#MAKE A LIST OF UMIs FROM THE UMI FILE
LIB_UMI_LIST = []
for line in LIB_UMI_LIST_FILE:
    line = line.strip().split()
    LIB_UMI_LIST.append(line[0])

#DEDUPLICATE
DUPLICATE_REF_DICT_FORWARD = {}
DUPLICATE_REF_DICT_REVERSE = {}
CHROM = 1
for line in INPUT_FILE:
    if line[0] == "@":
        OUTPUT_FILE.write(line)
    else:
        line = line.strip().split()
        umi = re.findall(":[A,G,C,T,N]{8}", line[0])[0].strip(":")
        if umi not in LIB_UMI_LIST:                                      #UNKNOWN UMI DON'T WRITE OUT
            continue
        else:
            flag = int(line[1])
            lposition = int(line[3])
            cigar = line[5]
            five_prime_pos, forward_mapped = cigar_parser(flag, cigar, lposition)
            if line[2] != CHROM:                                        #FLUSH DICTIONARY WHEN READS FROM A DIFFERENT CHROMOSOME IS ENCOUNTERED
                DUPLICATE_REF_DICT_FORWARD = {}
                DUPLICATE_REF_DICT_REVERSE = {}
            CHROM = line[2]
            ref_key = (CHROM, five_prime_pos)
            if forward_mapped == True:
                if ref_key not in DUPLICATE_REF_DICT_FORWARD:           #NOT A PCR DUPLICATE
                  DUPLICATE_REF_DICT_FORWARD[ref_key] = []
                  DUPLICATE_REF_DICT_FORWARD[ref_key].append(umi)
                  OUTPUT_FILE.write('\t'.join(line) + "\n")
                else:
                  if umi in DUPLICATE_REF_DICT_FORWARD[ref_key]:        #PCR DUPLICATE DON"T WRITE OUT
                    continue
                  else:
                    DUPLICATE_REF_DICT_FORWARD[ref_key].append(umi)     #NOT A PCR DUPLICATE
                    OUTPUT_FILE.write('\t'.join(line) + "\n")
            if forward_mapped == False:
                if ref_key not in DUPLICATE_REF_DICT_REVERSE:           #NOT A PCR DUPLICATE
                  DUPLICATE_REF_DICT_REVERSE[ref_key] = []
                  DUPLICATE_REF_DICT_REVERSE[ref_key].append(umi)
                  OUTPUT_FILE.write('\t'.join(line) + "\n")
                else:
                  if umi in DUPLICATE_REF_DICT_REVERSE[ref_key]:        #PCR DUPLICATE DON"T WRITE OUT
                    continue
                  else:
                    DUPLICATE_REF_DICT_REVERSE[ref_key].append(umi)
                    OUTPUT_FILE.write('\t'.join(line) + "\n")           #NOT A PCR DUPLICATE
            if forward_mapped == None:
               continue                                                 #UNMAPPED DON'T WRITE OUT
INPUT_FILE.close()
OUTPUT_FILE.close()
