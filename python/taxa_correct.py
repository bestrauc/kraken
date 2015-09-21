#!/usr/bin/python

import sys 

file_name = sys.argv[1]
tax_ids = []
tax_dict = {}
for arg in sys.argv[2:]:
    tax_ids.append(arg.split(','))

print tax_ids

handler = open(file_name, 'r')
for line in handler:
    count,tax = line.split()
    tax_dict[tax] = float(count)

classified=[]
s = sum(tax_dict.values())
classified.append(tax_dict['0'])

for tax in tax_ids:
    ts = 0
    for t in tax:
        ts = ts + tax_dict[t]
    classified.append(ts)

classified.insert(1,s - sum(classified))

for i in range(len(classified)):
    classified[i] = classified[i] / s

print classified
print sum(classified)
handler.close()
