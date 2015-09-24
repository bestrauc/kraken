#!/usr/bin/python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

if len(sys.argv) < 3:
    print "Usage: " + sys.argv[0] + " FILE LENGTHS"

handle = open(sys.argv[1], "rU")

for l in sys.argv[2:]:
    print "Trimming length %s." % l
    write_handle = open("%s_y.fasta" % l, "w")
    for record in SeqIO.parse(handle, "fastq"):
        new_rec = SeqRecord(record.seq[:int(l)], record.id, '', '')
	SeqIO.write(new_rec, write_handle, "fasta")
    handle.seek(0)
    write_handle.close()

handle.close()
#write_handle.close()



