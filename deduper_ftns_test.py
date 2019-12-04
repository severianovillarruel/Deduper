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

args = get_args()
INPUT_FILE = open(args.input_file, "r")
OUTPUT_FILE = open("dedupped.sam","w")
LIB_UMI_LIST_FILE = open("STL96.txt", "r")

######################################################################################

def flag_parser(flag):
  """Returns True if mapped 5' -> 3'. Returns False if mapped 3' -> 5'"""
  if((int(flag) & 4) != 4): #MAPPED
    if((int(flag) & 16) != 16): #5' -> 3' NORMAL (FORWARD) MAPPED
      return True
    else:
      return False #3' -> 5' REVERSE MAPPED
  else:
      return None  #NOT MAPPED (don't write out)

def cigar_parser(flag, cigar, lposition):
  """use the bitwise flag to see the direction the read mapped.
     Find the five prime position that the read mapped to.
     Find the number that is associated with a front end soft clipping if it exists.
     Add soft clipping at the beinning of the read to the position value to get the FIVE_PRIME_POS"""
  forward_map = flag_parser(flag)
  if forward_map == True:
      parsed_cigar = re.findall("^[0-9]+[S]", cigar)
      if parsed_cigar != []:
          parsed_cigar = int(parsed_cigar[0].strip("S"))
          five_prime_pos = (lposition - parsed_cigar)
      else:
          five_prime_pos = lposition
      return five_prime_pos, forward_map
  if forward_map == False:
      parsed_cigar = re.findall("[0-9]+", cigar)
      total_len = 0
      for num in parsed_cigar:
          num = int(num)
          total_len += num
      front_end_S_count = re.findall("^[0-9]+[S]", cigar)
      if front_end_S_count != []:
          front_end_S_count = int(front_end_S_count[0].strip("S"))
      else:
          front_end_S_count = 0
      insertion_count = re.findall("[0-9]+[I]", cigar)
      if insertion_count != []:
          insertion_count = int(insertion_count[0].strip("I"))
      else:
          insertion_count = 0
      five_prime_pos = total_len - front_end_S_count - insertion_count + lposition
      return five_prime_pos, forward_map
      #add everything except the soft clipping (when at the beginning of the CIGAR string) and the I values to the position to get the FIVE_PRIME_POS

######################################################################################

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
# print(DUPLICATE_REF_DICT_FORWARD)
print(DUPLICATE_REF_DICT_REVERSE)

#10S40M1I30M
INPUT_FILE.close()
OUTPUT_FILE.close()
