#!/usr/bin/python

import sys
import os 
import argparse
import struct
import subprocess
import time
import shutil
from collections import OrderedDict

# Dictionary of ranks and orders as they appear in nodes.dmp
order_dict = {'species':3, 'genus':2, 'family':1}

DATABASE_PATH = os.path.expanduser("~/kraken/kraken_db/database.kdb")
INDEX_PATH    = os.path.expanduser("~/kraken/kraken_db/database.idx")
TAXONOMY_PATH = os.path.expanduser("~/kraken/kraken_db/taxonomy/nodes.dmp")

#DATABASE_PATH = os.path.expanduser("~/kraken/minikraken_20140330/database.kdb")
#INDEX_PATH    = os.path.expanduser("~/kraken/minikraken_20140330/database.idx")
#TAXONOMY_PATH = os.path.expanduser("~/kraken/minikraken_20140330/taxonomy/nodes.dmp")

def tiles():
    """ generator for the tile numbers """
    for i1 in range(2):
        for i2 in range(3):
            for i3 in range(16):
                yield 1000*(i1+1) + 100*(i2+1) + i3+1

class KrakenFilter:

    # 'internal' functions not to be used from outside the class

    # Return values (X,Y) from X_Y strings (X,Y integers)
    def _splitReadID(self, ID):
        return [int(x) for x in ID.split('_')]

    def _reopenfilters(self):
        for tile, f in self.filters.iteritems():
            f = open(f.name, 'r+b')
            self.filters[tile] = f

    def _initializeFilters(self, path):
        self.filters = {}
        tile_offset = 0

        # read the tile files to determine number of reads
        for tile in tiles():
            tile_name = "s_1_"+str(tile)
            tile_handle = open(path+"/L001/C1.1/"+tile_name+".bcl", 'rb')
            N = struct.unpack('<i', tile_handle.read(4))[0]
            tile_handle.close()
 
            # Format: All 0, Version number (3), num. of clusters
            filter_header = struct.pack("<iii", 0x00000000, 0x00000003, N)
   
            # initialize filters with header and empty bytes        
            self.filters[tile] = open(path+"/L001/"+tile_name+".filter2",'wb') 
            self.filters[tile].write(filter_header)
            self.filters[tile].write(bytearray([0x01]*N))
            self.filters[tile].close()

            self.tile_index[tile] = tile_offset
            tile_offset += N
            print tile_offset

    def __init__(self, args, in_memory=True):
        self.tile_index = {}
        self._initializeFilters(args.BaseCalls)

        self.taxDB, self.taxTree = readTaxDB(args.taxonomy)
        self.specificity = args.specificity
        self.log_output_path = args.output
        self.log_output = open(self.log_output_path, 'w')
        self.classified_counter = 0

        ##self.reads = OrderedDict() if in_memory else None
        self.reads = [] if in_memory else None
       
        self.first_run = True
 
    # 'public' functions      

    # check if the taxID already satisfies our required specificity
    # or at least has an ancestor in the taxonomy with the specifity  
    def checkClassification(self, taxID):
        parent = taxID
        while (parent > 1) and (self.taxDB[parent] < self.specificity):
            parent = self.taxTree[parent]
        
        return (parent > 1)
    
#        if taxID in self.taxDB:
#            return self.taxDB[taxID] >= self.specificity
#        else: # Else case should only happen in case of a 0
#            return False;
    
    def filter_output(self, output_file, flush):      
        zero_byte = bytearray([0x00])
        self._reopenfilters()

        self.classified_counter = 0
        with open(output_file, 'r') as f:
            for line in f:
                fields = line.split('\t')
                tile, read = self._splitReadID(fields[1])
                taxID = int(fields[2])
           
                # add the reads to the list at the first run 
                if self.first_run and self.reads != None:
                    self.reads.append(line)
            
                if self.checkClassification(taxID):
                    self.filters[tile].seek(3*4 + read)
                    self.filters[tile].write(zero_byte)
                    self.classified_counter += 1

                    #self.reads[(tile,read)] = line
                    #self.log_output.write(line)
                
                # update the read with the most recent classification
                self.reads[self.tile_index[tile] + read] = line

        map(lambda h: h.close(), self.filters.values())
        #shutil.copyfile(self.log_output_path, self.log_output_path + str(l))
        
        self.first_run = False

        if flush:
            for line in self.reads:
                self.log_output.write(line)

    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("BaseCalls", help="The BaseCalls folder containing the lanes to be analyzed.")
    parser.add_argument("-o", "--output", help="Output file name for the final classification."
                        " (default: current directory)", default=os.getcwd()+'summary.class')
    parser.add_argument("--classifier", help="The path to the kraken classification binary.", 
                        metavar="path", required=True)
    parser.add_argument("-n", "--taxonomy", help="The path to the taxonomic tree (nodes.dmp) file.", 
                        metavar="path", required=True) 
    parser.add_argument("--specificity", help="The classification specifity at which we stop"
                        " classifying a read. (3-Species, 2-Genus, 1-Family. Default: 3)", default=2)
    parser.add_argument("-l", "--length", help="Length of the reads (number of cycle folders to be read)",
                        type=int, metavar="N", required=True)

    args = parser.parse_args()
        
    return args


def readTaxDB(tax_DB_file):
    taxDB = {}
    taxTree = {}
    with open(tax_DB_file, 'r') as f:
        for line in f:
            fields = line.split('\t|\t')
            tax_id = int(fields[0])
            tax_up = int(fields[1])
            rank = fields[2]
            taxTree[tax_id] = tax_up
            if rank in order_dict:
                taxDB[tax_id] = order_dict[rank]
            else:
                taxDB[tax_id] = 0

    return (taxDB,taxTree)

def main(argv):
    args = parse_args()

    root = os.path.abspath(args.BaseCalls)

    # Preliminary scan of the BaseCalls directory. 

    filter_files = {}
    krakenFilter = KrakenFilter(args)
   
    for lane_dir in os.listdir(root):
        lane_root  = os.path.join(root, lane_dir)
        lane_files = os.listdir(lane_root)
    
        filter_files[lane_dir] = [f for f in lane_files 
                        if os.path.splitext(f)[1] == ".filter"]

        # save at which cycle we start a new scan. 
        # Wait 31 cycles initially to get a 31-mer (Kraken default k-mer length)
        next_scan = 32
        # We wait until cycle N+1 starts so cycle N is fully written 
        next_scan_check_dir = os.path.join(lane_root, "C" + str(next_scan) + ".1")
         
        # while True because we want to do a do-while-loop
        while True:
            # monitor the directory until enough cycles for a new scan have been generated
            while not os.path.exists(next_scan_check_dir):
                print "Waiting"    
                time.sleep(5)

            print next_scan_check_dir, args.BaseCalls
            print [args.classifier,  "-d "+DATABASE_PATH, 
                   "-i "+INDEX_PATH, "-n "+args.taxonomy, 
                   "-t 4", "-b", "-M", "-l "+str(next_scan),
                   "-o tmp.class", args.BaseCalls]
            
            # start Kraken and analyze the directories up to next_scan
            # TODO: 'shell=True' security (shell injection) risk relevant here?
            subprocess.call(args.classifier + " -d "+DATABASE_PATH + 
                             " -i "+INDEX_PATH + " -n "+args.taxonomy + 
                             " -t 4" + " -b" + " -M" + " -l "+str(next_scan) +
                             " -o tmp.class " + args.BaseCalls, shell=True)

            # after Kraken finished runnning, scan the output to filter
            # reads that are classfied with sufficient specificity 

            krakenFilter.filter_output("tmp.class", flush=(next_scan==args.length));
            print krakenFilter.classified_counter
    
            if next_scan == args.length:
                break

            # wait until another 10 bases have been read before next scan
            next_scan = min(next_scan + 10, args.length)  
            next_scan_check_dir = os.path.join(lane_root, "C" + str(next_scan) + ".1")

    print filter_files
    
    print root
     

if __name__=="__main__":
    main(sys.argv)

