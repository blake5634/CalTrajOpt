#!/usr/bin/python3


#
#       verify how many unique paths are in a df
#

import numpy as np
import os
import c2to as cto
import sys
import brl_data.brl_data as bd

def main(args):

    #####################################################################
    #
    # ask user for a hint so they don't have to enter a long filename
    #
    datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    files = os.listdir(datadir)
    csvfiles = []
    for f in files:
        if '.csv' in f:
            csvfiles.append(f)

    hashstr = input('Enter first 4 characters of the hash from the filename:')
    if len(hashstr) < 4:
        print(' you like to live dangerously! (might match multiple files)')
    dfname = None
    for f in csvfiles:
        if hashstr in f:
            dfname = f
    if not dfname:
        bd.brl_error('Somethings wrong with hash string (not found)')

    #
    #  OK - now lets create a datafile and open it for reading
    #

    print('opening ', dfname)
    df = bd.datafile('', '','')  # open it with blank title info
    df.set_folders(datadir,'')        # set these to wherever you want to open datafiles
    df.open('r',tname=datadir+'/'+dfname)
    df.metadata.polish()  # convert metadata from strings to useful types
    print(df.metadata)

    cmin = df.metadata.d['Min Cost']
    cmax = df.metadata.d['Max Cost']

    x = input('pause ...')
    n=0
    pset = set()
    for row in df.reader:  # this is set up by the line df.open('r',tname=dfname)
        key = ''
        #row = df.type_row(row) # convert types
        for i in range(len(row)-1): #not including cost (last element)
            key+='{:3d}'.format(int(row[i])) # save time w.r.t type_row()
        print(' Adding:',key)
        pset.add(key)
        n+=1

    print('Unique paths: {:10}   rows:  {:10}'.format(len(pset),n))
    print('\n\n')




if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
