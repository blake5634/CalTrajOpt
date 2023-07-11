import c2to as cto
import brl_data.brl_data as bd
import sys
import os
import subprocess

# get list of files
datadir = ''
datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'

files = os.listdir(datadir)
mdflist = []
for f in files:
   if '_meta.json' in f:
      mdflist.append(str(f))

found = False
for f in mdflist:
   if sys.argv[1] in f:
      found = True
      subprocess.run(['kate', '--new' , datadir+'/'+f])

if not found:
   print('I found no metadata file with: ', sys.argv[1])
