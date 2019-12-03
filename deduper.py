#!/usr/bin/env python3
import re
import argparse
from deduper_ftns import *

def get_args():
    parser = argparse.ArgumentParser(description="dedupper")
    parser.add_argument("-file", "--input_file", help="Input file for deduplication", required=True, type = str)
    return parser.parse_args()

args = get_args()
INPUT_FILE = open(args.input_file, "r")
OUTPUT_FILE = open("dedupped.sam","w")
LIB_UMI_LIST_FILE = open("STL96.txt", "r")

#MAKE A LIST OF UMIs FROM THE LIBRARY
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
        if umi not in LIB_UMI_LIST:
            continue #don't write out
        else:
            flag = int(line[1])
            lposition = int(line[3])
            cigar = line[5]
            five_prime_pos, forward_mapped = cigar_parser(flag, cigar, lposition)
            if line[2] != CHROM:
                DUPLICATE_REF_DICT_FORWARD = {}
                DUPLICATE_REF_DICT_REVERSE = {}
            CHROM = line[2]
            ref_key = (CHROM, five_prime_pos)
            if forward_mapped == True:
                if ref_key not in DUPLICATE_REF_DICT_FORWARD:           #could modify to take the record with the highest qscore MAKE DIFFERENT DICTIONARY FOR FORWARD AND REVERSE
                  DUPLICATE_REF_DICT_FORWARD[ref_key] = []
                  DUPLICATE_REF_DICT_FORWARD[ref_key].append(umi)
                  OUTPUT_FILE.write('\t'.join(line) + "\n")
                else:
                  if umi in DUPLICATE_REF_DICT_FORWARD[ref_key]:
                    continue        #don't print out
                  else:
                    OUTPUT_FILE.write('\t'.join(line) + "\n")
            if forward_mapped == False:
                if ref_key not in DUPLICATE_REF_DICT_REVERSE:           #could modify to take the record with the highest qscore MAKE DIFFERENT DICTIONARY FOR FORWARD AND REVERSE
                  DUPLICATE_REF_DICT_REVERSE[ref_key] = []
                  DUPLICATE_REF_DICT_REVERSE[ref_key].append(umi)
                  OUTPUT_FILE.write('\t'.join(line) + "\n")
                else:
                  if umi in DUPLICATE_REF_DICT_REVERSE[ref_key]:
                    continue        #don't print out
                  else:
                    OUTPUT_FILE.write('\t'.join(line) + "\n")
            if forward_mapped == None:
               continue     #don't write out

INPUT_FILE.close()
OUTPUT_FILE.close()
