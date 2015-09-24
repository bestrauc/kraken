#!/usr/bin/python

import sys
import os 
import struct
import binascii
import argparse
import random

#import Bio 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

DNAdict = {'A':0, 'C':1, 'G':2, 'T':3}

def tiles():
    """ generator for the tile numbers """
    for i1 in range(2):
        for i2 in range(3):
            for i3 in range(16):
                yield 1000*(i1+1) + 100*(i2+1) + i3+1

def open_file_handles(n, tile, base_path):
    handles = []
    for i in range(n):
        h = open(base_path + "/BaseCalls/L001/C%i.1/s_1_%i.bcl"%((i+1), tile), "w+b")
        h.write(struct.pack("<i",0)) # write placeholder 32bits
        handles.append(h)
    
    return handles

def create_directories(n, base_path):
    for i in range(n):
        os.makedirs(base_path + "/BaseCalls/L001/C%i.1"%(i+1))

# Wrapper we use to modify the SeqIO.parse generator so that it
# returns a fixed number (random_num) of random reads apart from
# the original number (read_num) of fastq reads.
def SeqIO_addRandom_fastq(SeqIO_fastq_parser, read_num, random_num, read_length):
    def RandomDNA():
        return ''.join(random.choice('ACGT') for _ in range(read_length))

    total = read_num + random_num
    for i in range(total):
        # probability of yielding a true read
        p = float(read_num)/float(read_num + random_num)

        # yield true read from file
        if random.random() < p:
            read_num = read_num - 1
            yield (1, SeqIO_fastq_parser.next())
        # yield simulated read
        else:    
            random_num = random_num - 1
            qual_annotation = {"phred_quality": [0]*read_length}
            random_record = SeqRecord(Seq(RandomDNA(), generic_dna),
                                      id="random", 
                                      letter_annotations=qual_annotation)
            
            yield (0, random_record)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="The fastq-file to generate the tile files from.")
    parser.add_argument("-o", "--output", help="Output directory for BaseCalls folder."
                        " (default: current directory)", default=os.getcwd())
    parser.add_argument("--id", help="Write sequence IDs into text files.", 
                        action="store_true", default=True)
    parser.add_argument("--tiles", help="Number of tiles to generate. (Default: 96)", 
                        type=int, choices=range(1,97), metavar="[1-96]", default=96)
    parser.add_argument("--filtered", help="Fraction of reads that did not pass the filter."
                        " (Default: 0)", type=float, default=0)
    parser.add_argument("--add_noise", help="If true, filtered reads are randomly generated.",
                        action="store_true", default=False)
    args = parser.parse_args()
        
    return args

def main(argv):
    args = parse_args() 

    # Count the number of sequences to ensure equal distribution among tiles.
    num_lines = sum(1 for line in open(args.file))
    assert(num_lines % 4 == 0) # maybe use if-clause instead of assert

    read_count = num_lines / 4
    random_size = int((args.filtered*read_count)/(1-args.filtered))

    count = read_count
    # add the number of false noise reads to the count, if we add those
    if args.add_noise:
        count = read_count + random_size

    # Read the length of the fastq reads (we only read one record)
    handle = open(args.file, "rU")
    for record in SeqIO.parse(handle,"fastq"):
        length = len(record.seq)
        break 
    handle.close()

    # Set tile number generator and reads per tile.
    tile_size = count / args.tiles
    tile = enumerate(tiles())

    # Format: All 0, Version number (3), All 0 placeholder for num. of clusters
    filter_header = struct.pack("<iii", 0x00000000, 0x00000003, 0x00000000)

    # Create the BaseCalls/LYY/CX.1 subdirectory structure
    create_directories(length, args.output)
    t_i, active_tile = tile.next()
    handles = open_file_handles(length, active_tile, args.output)

    read_records = 0
    input_handle = open(args.file, "r")
    filter_handle = open(args.output + "/BaseCalls/L001/s_1_%i.filter"%active_tile, "w")
    filter_handle.write(filter_header)
    
    if args.id:
        name_handle = open(args.output + "/BaseCalls/L001/s_1_%i.names"%active_tile, "w")
 
    # Output progress status.
    sys.stdout.write("Writing tile: %d" % (active_tile)) 
    sys.stdout.flush()

    SeqIO_iter = SeqIO.parse(input_handle, "fastq")   

    # wrap the SeqIO parser to generate random reads (none if no noise should be added)
    if args.add_noise:
        SeqIO_iter = SeqIO_addRandom_fastq(SeqIO_iter, read_count, random_size, length)
    else:
        SeqIO_iter = SeqIO_addRandom_fastq(SeqIO_iter, read_count, 0, length)

    for filtered, record in SeqIO_iter:
        if (read_records == tile_size and (t_i+1) < args.tiles):
            # Write the first byte, now we know how many reads were written.
            map(lambda handle: handle.seek(0), handles)
            map(lambda handle: handle.write(struct.pack("<i", read_records)), handles)
            map(lambda handle: handle.close(), handles)

            # Update read count
            filter_handle.seek(8)
            filter_handle.write(struct.pack("<i", read_records))
            filter_handle.close()

            # Prepare next iteration.            
            read_records = 0
            t_i, active_tile = tile.next()
            
            handles = open_file_handles(length, active_tile, args.output)
            filter_handle = open(args.output + "/BaseCalls/L001/s_1_%i.filter"%active_tile, "w")
            filter_handle.write(filter_header)           
 
            if args.id:
                name_handle = open(args.output + "/BaseCalls/L001/s_1_%i.names"%active_tile, "w")
 
            
            # Output progress status.
            sys.stdout.write("\rWriting tile: %d" % (active_tile)) 
            sys.stdout.flush()
     
        qual = record.letter_annotations["phred_quality"]

        # Write sequence to tile files.
        for i,c in enumerate(record.seq):
            if (c=='N'):
                data = struct.pack("B", 0)
            else:
                data = struct.pack("B", (qual[i] << 2) | DNAdict[c])

            handles[i].write(data)
   
        if args.add_noise: 
            filter_handle.write(struct.pack("B", filtered))
        else:
            filter_handle.write(struct.pack("B", random.random() >= args.filtered))

        if args.id:
            name_handle.write(record.id + "\n")

        read_records += 1

    # Write the first sequence number byte for the last tile.
    map(lambda handle: handle.seek(0), handles)
    map(lambda handle: handle.write(struct.pack("<i",read_records)), handles)
    map(lambda handle: handle.close(), handles)

    # Update read count
    filter_handle.seek(8)
    filter_handle.write(struct.pack("<i", read_records))
    filter_handle.close()

if __name__=="__main__":
    main(sys.argv)

# TODO: Add .stats etc. files to output.
