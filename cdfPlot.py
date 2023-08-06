import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import sys
import os

import brl_data.brl_data as bd


def main(args):
    ###################################################################
    #
    #ask user for a hint so they don't have to enter a long filename
    #
    datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'

    #  all the dirs we might find data files
    dirs = ['/home/blake/Sync/Research/CalTrajOpt_RESULTS',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data/Gold'
            ]

    fnames = []
    files = []
    for d in dirs:
        fnames = os.listdir(d)
        for n in fnames:
            files.append([d,str(n)]) # dir, name

    csvfiles = []
    for f in files: # [0] = dir [1] = name
        if '.csv' in f[1]:
            csvfiles.append(f)

    #hashstr = input('Enter first 4 characters of the hash from the filename:')
    hashstr = args[1]
    if len(hashstr) < 4:
        print(' you like to live dangerously! (might match multiple files)')
    dfname = None
    tieflag = 'ties_info'
    for f in csvfiles:
        if hashstr in f[1] and tieflag not in f[1]:
            dfname = f[1]
            dfdir  = f[0]
    if not dfname:
        bd.brl_error('Somethings wrong with hash string (not found)')


    print('opening ', dfdir + '/'+ dfname)
    df = bd.datafile('', '','')  # open it with blank title info
    df.set_folders(dfdir,'')        # set these to wherever you want to open datafiles
    df.open('r',tname=dfdir+'/'+dfname)
    df.metadata.polish()  # convert metadata from strings to useful types


    databinvals = df.metadata.d['CostHistogram_values']
    databinlevels = df.metadata.d['CostHistogram_levels']
    mu = df.metadata.d['CostMean']
    sd = df.metadata.d['CostStDev']

    totalsum = sum(databinvals)
    prop0 = databinlevels[0]
    prop1 = databinlevels[-1]
    dom = prop1-prop0

    # normalize:
    normdblevels = []
    normdbvals = []
    ncdf = []
    for dl in databinlevels:
        ndl = (dl-mu)/sd
        normdblevels.append(ndl)
        ncdf.append(stats.norm.cdf(ndl))
    for dv in databinvals:
        normdbvals.append(dv/totalsum)

    #print('dbv:')
    #print(normdbinvals)
    #print('dblev:')
    #print(databinlevels)
    cumdist = [0.0]
    tot = 0.0
    for dv in normdbvals:
        tot += dv
        cumdist.append(tot)


    my_dpi = 200
    w = 1250
    h = 1250
    plt.figure(figsize=(w/my_dpi, h/my_dpi), dpi=my_dpi)
    plt.plot(normdblevels, cumdist)
    plt.plot(normdblevels, ncdf) #, databinlevels,databinvals)
    plt.show()



if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
