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
    # ask user for a hint (argv[1]) so they don't have to enter a long filename
    #
    datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    dirlist = [datadir, '/home/blake/Ptmp/CalTrajOpt']

    mf = bd.finder()
    mf.set_dirs(dirlist)

    hashstr = args[1]
    res = mf.findh([hashstr, '.csv']) # find files with hash and are .csv ext.

    if len(res) < 1:
        brl_error(' file not found. hash: '+hashstr)
    if len(res) > 1:
        brl_error('   too many results: '+hashstr)
    datadir = res[0][0]
    dfname  = res[0][1]

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

    sum1 = 0.0  # sum of c
    sum2 = 0.0  # sum of c^2
    n = 0
    for row in df.reader:  # this is set up by the line df.open('r',tname=dfname)
        df.print_row(row)  # automatically formats the row according to types in metadata
        c = row[-1]
        for i,bl in enumerate(binlims):
            if c < bl:
                bins[i] += 1
                break
        #accumulate sums and count
        sum1 += c
        sum2 += c*c
        n += 1
    mean = sum1/n
    sd   = np.sqrt((n*sum2-sum1*sum1)/(n*(n-1)))  # sample sd

    print('Mean: {:5.2f}   sd:  {:5.2f}'.format(mean,sd))
    df.metadata.d['CostMean'] = mean
    df.metadata.d['CostStDev'] = sd

    print('Cost histogram: {:}'.format(df.metadata.d['CostType']))
    fullims = np.append(np.array([cmin]),binlims)

    print(' ... please see metadata: ', df.hashcode)
    #for i,b in enumerate(bins):
        #if i>0:
            #print('{:3d} {:5.2f} - {:5.2f} | {:}'.format(i,fullims[i-1],fullims[i],b[i-1]))

    df.metadata.d['CostHistogram_values'] = bins
    df.metadata.d['CostHistogram_levels'] = list(fullims)
    df.metadata.save()
    df.close()

    if False:  #CSV version
        print('\n\n    CSV version\n')
        for i,b in enumerate(bins):
            if i>0:
                print('{:5.2f}, {:}'.format(fullims[i-1],b[i-1]))
    print('\n\n')




if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
