#!/usr/bin/python3


#
#       Collect histogram of # paths vs cost
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

    print('')
    print('first rows of data are:')
    df.print_header()
    n=0

    cmin = df.metadata.d['Min Cost']
    cmax = df.metadata.d['Max Cost']
    NBINS = 20
    step = (cmax-cmin)/NBINS
    binlims = np.arange(cmin,cmax,step)
    print('Cmin: {:4.2f}  Cmax: {:4.2f}'.format(cmin,cmax))
    print('working with bin limits: ',binlims)
    bins = [0]*NBINS
    binlims = binlims[1:] # chop the first value and add last value
    binlims = np.append(binlims,cmax)
    print('working with bin limits: ',binlims)
    x = input('pause ...')
    for row in df.reader:  # this is set up by the line df.open('r',tname=dfname)
        df.print_row(row)  # automatically formats the row according to types in metadata
        c = row[9]
        for i,bl in enumerate(binlims):
            if c < bl:
                bins[i] += 1
                break

    print('Cost histogram: {:}'.format(df.metadata.d['CostType']))
    fullims = np.append(np.array([cmin]),binlims)
    for i,b in enumerate(bins):
        if i>0:
            print('{:3d} {:5.2f} - {:5.2f} | {:}'.format(i,fullims[i-1],fullims[i],b))
    #df.metadata.d['CostHistogram_values'] = bins
    #df.metadata.d['CostHistogram_levels'] = list(fullims)
    df.metadata.d['CostHistogram_values_levels'] = zip(bins,list(fullims))
    df.close()

    if True:  #CSV version
        print('\n\n    CSV version\n')
        for i,b in enumerate(bins):
            if i>0:
                print('{:5.2f}, {:}'.format(fullims[i-1],b))
    print('\n\n')




if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
