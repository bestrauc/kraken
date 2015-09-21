#!/usr/bin/python

import sys
import argparse

print sys.argv

parser = argparse.ArgumentParser(description='Check if the taxonomic classification changed ' 
                                             'between Kraken runs.')
parser.add_argument('taxonomy', help='Tax-IDs of the taxa we expect, ascending in specifity. '
                                     '(e.g. Family -> Genus -> Species)', nargs='+')
parser.add_argument('--files', help='Kraken output files to analyze.', required=True)
#parser.add_argument('--taxonomy', help='Tax-IDs of the taxa we expect, ascending in specifity. '
#                                       '(e.g. Family -> Genus -> Species)')

args = parser.parse_args()

_file = open(args.files,"r")
kraken_files = [ line.strip() for line in _file ]
_file.close()

print args

# create dictionary filled with tax ids.
read_dict = {}

taxonomy = args.taxonomy

for f in kraken_files:
  handle=open(f,'r')  
  for line in handle:
    tmp = line.split()
    if tmp[1] not in read_dict:
        read_dict[tmp[1]] = []

    read_dict[tmp[1]].append(tmp[2])

# detect changes in tax ids
changes_dict = {}
for key in read_dict:
    changes_dict[key] = 0
    last_val = read_dict[key][0] 
    for val in read_dict[key][1:]:
        if (last_val in taxonomy and val not in taxonomy) or \
           (last_val in taxonomy and val in taxonomy \
            and taxonomy.index(last_val) > taxonomy.index(val)):
            changes_dict[key] += 1
        last_val = val

    if changes_dict[key] > 0:
        print read_dict[key]

vals = changes_dict.values()
print float(sum(vals))/len(vals)
