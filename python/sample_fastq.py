#!/usr/bin/python

import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("Ref", help="File containing paths to reference files.")
    parser.add_argument("-l", "--length", help="Length of the sampled reads.", 
                        default=100)
    parser.add_argument("-o", "--output", help="Output file name for the .fastq file.",
                        required=True)

    args = parser.parse_args()
        
    return args

def main(argv):
    args = parse_args()

    files = [line.strip().split() for line in open(args.Ref).readlines()]

    outfile = open(args.output, "w")   
 
    for f in files[:10]:
        #TODO: Does not work if reference file folder is read only
        command = "mason_simulator -ir %s --illumina-read-length %s -n 100000 -o %s" \
                    %(f[0], args.length, "tmp.fastq")
        print command
        os.system(command)

        with open("tmp.fastq") as tmp_in:
            for line in tmp_in:
                outfile.write(line)

    os.remove("tmp.fastq")

if __name__=="__main__":
    main(sys.argv)
