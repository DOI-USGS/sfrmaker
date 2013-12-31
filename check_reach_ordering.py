#  check_reach_ordering.py
#  short script to read in reach_ordering.txt
#  and check to see that every downstream hydrosequence number
#  that is not zero has a corresponding hydrosequence entry
#
#  hydrosequence points to the next downstream hydrosequence
#  the next downstream one should exist in the data file....
#
import os
import re
import sys

dnhydroseq=dict()
hydroseq=dict()

ORDER=open('reach_ordering.txt','r')

header=ORDER.readline()

for line in ORDER:
    vals=re.split(',',line)
    hydroseq[int(vals[2])]=1
    dnhydroseq[int(vals[4])]=1

ORDER.close()

for dwn in dnhydroseq.iterkeys():
    if not dwn in hydroseq and dwn>0:
        print '%d not in hydroseq' % dwn
        
