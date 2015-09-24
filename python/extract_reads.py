# test script to extract reads from one tile only and only the first 100 bp
import struct

def bcl(lane,cycle,tile=1101):
    return 'L00%i/C%i.1/s_%i_%04i.bcl'%(lane,cycle,lane,tile)

def tiles():
    """ generator for the tile numbers """
    for i1 in range(2):
        for i2 in range(3):
            for i3 in range(16):
                yield 1000*(i1+1) + 100*(i2+1) + i3+1

def tile_id(tile):
    """ unique id between 0 and 96 for each tile """
    t = int(tile)
    i1 = t/1000
    i2 = (t-1000*i1)/100
    i3 = t - 1000*i1 - 100*i2
    return (i1-1)*48 + (i2-1)*16 + i3-1



for lane in range(1,9):
    print 'creating data structures'
    # data structures for all reads on this lane
    reads = [0]*96
    # get the number of reads on each tile from the first cycle
    for tile in tiles():
        fn = bcl(lane,1,tile)
        f = open(fn,'rb')
        # number of reads in this bcl file
        size = struct.unpack('<i',f.read(4))[0]
        reads[tile_id(tile)] = ['']*size
        f.close()
   
    # now read all reads for the lane
    for cl in range(1,2):
        print 'reading cycle',cl
        for tile in tiles():
            tid = tile_id(tile)
            fn = bcl(lane,cl,tile)
            f = open(fn,'rb')
            # number of reads in this bcl file
            size = struct.unpack('<i',f.read(4))[0]
            # read bunch
            for i in range(size):
                raw = ord(f.read(1))
                base = int(raw & 0b11 )
                """if base == 0:
                    reads[tid][i] += 'A'
                elif base == 1:
                    reads[tid][i] += 'C'
                elif base == 2:
                    reads[tid][i] += 'G'
                elif base == 3:
                    reads[tid][i] += 'T'"""
            f.close()
            
    print len(reads)
    break

import sys
sys.exit(1)

for cl in range(1,101):
    # open bcl file
    fn = bcl(1,cl,1101)
    f = open(fn,'rb')
    # number of reads in this bcl file
    size = struct.unpack('<i',f.read(4))
    # read bunch
    for i in range(size):
        raw = ord(f.read(1))
        base = int(raw & 0b11 )
        if base == 0:
            reads[i] += 'A'
        elif base == 1:
            reads[i] += 'C'
        elif base == 2:
            reads[i] += 'G'
        elif base == 3:
            reads[i] += 'T'
        else:
            raise Exception('Error: not recognized letter: %i'%base)

print reads[0]


f_out = open('test.fasta','w')
for i,r in enumerate(reads):
    f_out.write('>seq_%i\n'%i)
    f_out.write(r+'\n')
f_out.close()
