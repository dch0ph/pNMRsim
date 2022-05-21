#!/usr/bin/python
import sys
from Numeric import *

nargs=len(sys.argv)

if nargs<2:
    print 'Syntax: makeexchange <input file>'
    sys.exit()

f=file(sys.argv[1])

matsize=0
lineno=1
rowcount=0
for line in f.xreadlines():
    if not line.startswith('#'):        
        els=line.split()
        nels=len(els)
        if matsize==0:
            matsize=nels
            exchmat=zeros([matsize,matsize],Float)
        elif matsize != nels:
            print 'Number of columns changed from', matsize, ' to', nels, ' at line', lineno, '!'
            sys.exit()
        whichrow = rowcount % matsize
        for i in range(matsize):
            exchmat[whichrow][i]=exchmat[whichrow][i]+float(els[i])
        rowcount=rowcount+1
    lineno=lineno+1

f.close()
if rowcount % matsize:
    print 'Warning: number of data rows ', rowcount, ' is not divisible by size of exchange matrix', matsize

for r in range(matsize):
    for c in range(matsize):
        sys.stdout.write(str(exchmat[r][c])+'\t')
    sys.stdout.write('\n')


