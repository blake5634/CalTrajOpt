#
import os
import sys
from pathlib import Path
import brl_data.brl_data as bd

dir = '/home/blake/Ptmp/CalTrajOpt/testing/'
keeperfiles = []
targetfiles = []

# 4 target files and 4 keeper files (stored URIs)
for i in range(4):
   targetfiles.append(f'target{bd.brl_id(8)}.txt')

f = open(f'{dir}keepers.txt','r')
for URI in f:
   keeperfiles.append(f'keeper{URI.strip()}.txt')

for f in ['folder1/','folder2/','folder3/']:
   # clean out old files
   f2 = dir+f
   filelist = os.listdir(f2)
   for fn in filelist:
      print(f'cleaning up {os.path.join(f2,fn)}')
      os.remove(os.path.join(f2,fn))
   # create new ones
   for k in keeperfiles:
      keeper = dir+f+k
      print('creating: ',keeper)
      Path(keeper).touch()
   for t in targetfiles:
      target = dir+f+t
      Path(target).touch()

