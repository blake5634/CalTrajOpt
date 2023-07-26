#!/usr/bin/python3

import c2to as cto
import brl_data.brl_data as bd
import sys
import os
import subprocess

#  all the dirs we might find data files
dirs = ['/home/blake/Sync/Research/CalTrajOpt_RESULTS',
        '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data',
        '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D-Round2',
        '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data/Gold'
        ]

datadir = 'Using list'
print('cmd line: ',sys.argv)
if len(sys.argv) < 2:
   bd.brl_error('To view metadata:\n  usage: >mdview FILE ... or ... > mdview dir FILE')
# get list of files
if len(sys.argv) > 2:
   datadir = sys.argv[1]
   targethash = sys.argv[2]
elif len(sys.argv) == 2:
   targethash = sys.argv[1]

print('datadir: ', datadir)
print('targethash: ', targethash)

fnames = []
files = []
for d in dirs:
   fnames += os.listdir(d)
   for n in fnames:
      files.append([d,str(n)]) # dir, name

mdflist = []
for f in files:
   if '_meta.json' in f[1] and targethash in f[1]:
      mdflist.append(f)

found = False
if len(mdflist) > 1:
   mdflist = [mdflist[0]]    # de dup (copies in other dirs)
for f in mdflist:
   if targethash in f[1]:
      found = True
      filepath = f[0]+'/'+f[1]
      subprocess.run( ['kate', '--new', filepath] )

if not found:
   print('I found no metadata file with: ', targethash ,'in:')
   for d in dirs:
      print('     ',d)
