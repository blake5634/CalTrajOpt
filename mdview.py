#!/usr/bin/python3

import c2to as cto
import brl_data.brl_data as bd
import sys
import os
import subprocess

print('cmd line: ',sys.argv)
if len(sys.argv) < 2:
   bd.brl_error('To view metadata:\n  usage: >mdview FILE ... or ... > mdview dir FILE')
# get list of files
if len(sys.argv) > 2:
   datadir = sys.argv[1]
   targethash = sys.argv[2]
elif len(sys.argv) == 2:
   targethash = sys.argv[1]
   datadir = ''
   datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'

print('datadir: ', datadir)
print('targethash: ', targethash)

files = os.listdir(datadir)
mdflist = []
for f in files:
   if '_meta.json' in f:
      mdflist.append(str(f))

found = False
for f in mdflist:
   if targethash in f:
      found = True
      subprocess.run(['kate', '--new' , datadir+'/'+f])

if not found:
   print('I found no metadata file with: ', datadir+'/ ... '+targethash)
