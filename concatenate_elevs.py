indat = []
for i in range(19):
    tmp = open('elev%d.dat' %(i+1),'r').readlines()
    indat.extend(tmp)
    print i+1
    
ofp = open('allbottoms.dat','w')
for line in indat:
    ofp.write(line.strip() + '\n')
ofp.close()