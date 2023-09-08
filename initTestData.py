#
import os
import shutil
import sys
from pathlib import Path
import brl_data.brl_data as bd

tdir = '/home/blake/Ptmp/CalTrajOpt/testing/'

logs = ['/home/blake/Ptmp/CalTrajOpt/testing/log1.txt',
        '/home/blake/Ptmp/CalTrajOpt/testing/log2.txt' ]
testfolders = ['folder1/','folder2/','folder3/']
trash = '/home/blake/Ptmp/CalTrajOpt/testing/Trash'

try:
   shutil.rmtree(tdir)
except:
   pass
os.mkdir(tdir)
for fld in testfolders:
   os.mkdir(tdir+fld)
os.mkdir(trash)

f = open(tdir+'keepers.txt', 'w')
for l in ['aaaaaaaa', '1234abcd', 'cc2cff5f', 'dead4568']:
   print(l,file=f)
f.close()

keeperfiles = []
targetfiles = []
keeperURIs  = []

# 4 target files and 4 keeper files (stored URIs)
for i in range(4):
   targetfiles.append(f'target{bd.brl_id(8)}')

f = open(f'{tdir}keepers.txt','r')
for URI in f:
   keeperURIs.append(URI.strip())
   keeperfiles.append(f'keeper{URI.strip()}')

for i,f in enumerate( testfolders ):
   # clean out old files
   f2 = tdir+f
   filelist = os.listdir(f2)
   for fn in filelist:
      print(f'cleaning up {os.path.join(f2,fn)}')
      os.remove(os.path.join(f2,fn))
   # create new ones
   for k in keeperfiles:
      keeper = tdir+f+k+f'_f{i:02}'+'.txt'
      print('creating: ',keeper)
      Path(keeper).touch()
   for t in targetfiles:
      target = tdir+f+t+f'_f{i:02}'+'.txt'
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
