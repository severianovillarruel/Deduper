#!/usr/bin/env python3
import re

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
