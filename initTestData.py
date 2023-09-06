#
import os
import sys
from pathlib import Path
import brl_data.brl_data as bd

dir = '/home/blake/Ptmp/CalTrajOpt/testing/'

logs = ['/home/blake/Ptmp/CalTrajOpt/testing/log1.txt',
        '/home/blake/Ptmp/CalTrajOpt/testing/log2.txt' ]

keeperfiles = []
targetfiles = []
keeperURIs  = []

# 4 target files and 4 keeper files (stored URIs)
for i in range(4):
   targetfiles.append(f'target{bd.brl_id(8)}')

f = open(f'{dir}keepers.txt','r')
for URI in f:
   keeperURIs.append(URI.strip())
   keeperfiles.append(f'keeper{URI.strip()}')

for i,f in enumerate( ['folder1/','folder2/','folder3/'] ):
   # clean out old files
   f2 = dir+f
   filelist = os.listdir(f2)
   for fn in filelist:
      print(f'cleaning up {os.path.join(f2,fn)}')
      os.remove(os.path.join(f2,fn))
   # create new ones
   for k in keeperfiles:
      keeper = dir+f+k+f'_f{i:02}'+'.txt'
      print('creating: ',keeper)
      Path(keeper).touch()
   for t in targetfiles:
      target = dir+f+t+f'_f{i:02}'+'.txt'
      print('creating: ',target)
      Path(target).touch()

#
#   generate log files with some lines (targets) which do not have URLs in the list
#          and also some log entries without any URLs at all (save these)
for lf in logs:
   n = len(keeperURIs)
   f = open(lf,'w')
   # generate a log with keepers and targets
   for i in range(n):
      le_k = f'Some text describing a file we want to keep {keeperURIs[i]}'
      le_t = f'Some text describing a file we want to ACTION {bd.brl_id(8)}'
      le_n =  'Some text describing without a URI at all (a note)'
      print(le_k,file=f)
      print(le_t,file=f)
      print(le_n,file=f)
   f.close()
