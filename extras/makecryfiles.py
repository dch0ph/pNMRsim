#!/usr/bin/python
import sys
#from Numeric import *

nargs=len(sys.argv)

if nargs<2:
    print 'Syntax: makeexchange <input file>'
    sys.exit()

fname=sys.argv[1]
f=file(fname)

matsize=0
lineno=1
rowcount=0
orientinfo=[]
masterstack=[]
donematrix=False
for line in f.xreadlines():
    if not line.startswith('#'):        
        els=map(float,line.split())
        nels=len(els)-4
        if matsize==0:
            matsize=nels
        elif matsize != nels:
            print 'Number of columns changed from', matsize, ' to', nels, ' at line', lineno, '!'
            sys.exit()
	if rowcount==0:
	    exchmat=[]
	    orient=els[matsize:matsize+3]
	    weight=els[matsize+3]
        exchmat.extend(els[:matsize])
        rowcount=rowcount+1
	donematrix= (rowcount==matsize)
	if donematrix:
		masterstack.append([orient,weight,exchmat])
		rowcount=0
    lineno=lineno+1

f.close()
if not donematrix:
    print 'Warning: incomplete or missing exchange matrix'

#print masterstack
for r in range(matsize):
    for c in range(matsize):
	ind=r*matsize+c
	newfname= fname + str(r+1) + ',' + str(c+1)
	f=open(newfname,'w')
	for curitem in masterstack:
		for angle in curitem[0]:
			f.write(str(angle)+'\t ')
		weight=curitem[1]
		el=curitem[2][ind]
		f.write(str(weight*el)+'\n')
	f.close()



