import glob
import os
process=7
shellloc = 'C:/cygwin64/home/YK/PyWeath/'
seriesname='test'
filelist = glob.glob(shellloc+seriesname+'*.sh')
for i in range(len(filelist)):
    os.remove(filelist[i])
wholework=shellloc+seriesname+'.sh'
line = '#!/bin/bash -e\n'\
       + 'a=1\n'\
       + 'while [ $a -lt ' +str(process+1) +' ]\n'\
       + 'do \n'\
       + '  echo $a\n'\
       + '  dos2unix test_${a}.sh\n'\
       + '  chmod u+x "'+seriesname+'_${a}.sh"\n'\
       + '  nohup ./'+seriesname+'_${a}.sh' \
       + '> out${a}.log < /dev/null &\n' \
       + '  a=`expr $a + 1`\n' \
       + 'done\n'
f=open(wholework, 'w')
f.write(line)
f.close()

cnt=0
for i in range(16):
    for j in range(11):
        
        runfile = shellloc + seriesname +'_'+str(cnt%process+1)+'.sh'
        f=open(runfile,'a')
        runplace = './run/run'+str(cnt+1)
        line = runplace \
                   +'\n'
        f.write(line)
        f.close()
        cnt+=1
