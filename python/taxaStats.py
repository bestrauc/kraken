#!/usr/bin/python

import sys
import argparse

order_dict = {'species': 3, 'genus': 2, 'family':1}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("kraken_ref", help="A kraken reference classification output file.")
    parser.add_argument("kraken_inc", help="An incremental kraken classification output file.")
    parser.add_argument("-t", "--taxonomy", help="Taxonomy file containing tax-id ranks.",
                        required=True)

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

def on_path(taxTree, taxid1, taxid2):
    node = taxid1    
    while (node!=taxid2 and node!=1):
        node = taxTree[node]

    return (node == taxid2)
    

def main():
    args = parse_args()   

    f1 = open(args.kraken_ref)
    f2 = open(args.kraken_inc)

    taxDB, taxTree = readTaxDB(args.taxonomy)

#    counter = {3:0, 2:0, 1:0, 0:0}
    equals_count = 0
    path_count = 0

    for line1 in f1:
        line2  = f2.readline()
        fields1 = line1.split('\t')
        fields2 = line2.split('\t')
        taxid1  = int(fields1[2])
        taxid2  = int(fields2[2])

        if taxid1 == taxid2:
            equals_count += 1

        # check whether the 2nd classification lies on the path to the root of the 1st
        # if the ancestor satisfies our classification requirement we count it
       
        if on_path(taxTree, taxid1, taxid2):
            path_count += 1
        else:
            print line1
            print line2
            print '======================================'

    print equals_count
    print path_count


if __name__=="__main__":
    main()
