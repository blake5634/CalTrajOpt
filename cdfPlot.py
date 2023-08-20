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
<<<<<<< HEAD
    goalhash =   args[1]
    myff = bd.finder()
    myff.set_dirs(dirs)
    fds = myff.findh([ goalhash ])

    csvfiles = []
    for f in fds: # [0] = dir [1] = name
        if '.csv' in f[1]:
            csvfiles.append(f)

    if len(csvfiles) < 1:
        bd.brl_error(f'Somethings wrong with hash string {goalhash} (not found)')
    else:
        tieflag = 'ties_info'
        for f in csvfiles:
            if goalhash in f[1] and tieflag not in f[1]:
                dfdir  = f[0]
                dfname = f[1]

    print('opening ', dfdir + '/'+ dfname)
    fullhash = bd.getHashFromFilename(dfname)

    df = bd.datafile('', '','')  # open it with blank title info
    df.set_folders(dfdir,'')        # set these to wherever you want to open datafiles
    df.open('r',tname=dfdir+'/'+dfname)
    df.metadata.polish()  # convert metadata from strings to useful types

    ########################################################################
    # convert the histogram to sample cdf
    databinvals = df.metadata.d['CostHistogram_values']
    databinlevels = df.metadata.d['CostHistogram_levels']
    mu = df.metadata.d['CostMean']
    sd = df.metadata.d['CostStDev']

    totalsum = sum(databinvals)
    prop0 = databinlevels[0]
    prop1 = databinlevels[-1]
    dom = prop1-prop0

    #############################
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

    cumdist = [0.0]
    tot = 0.0
    for dv in normdbvals:
        tot += dv
        cumdist.append(tot)

    ##########################################################################
    #  Try the ks-static
    #
    ncdata = []
    for row in df.reader:
        c = float(row[-1])
        ncdata.append((c-mu)/sd)

    k = stats.ks_1samp(ncdata, stats.norm.cdf)

    ChiSq = k.statistic
    pval = k.pvalue

    print('Koloogorov-Smirnov results: ChiSq: {:8.5f}  pval = {:8.5f}, n={:}'.format(ChiSq, pval, len(ncdata)))

    ########################
    #  Q-Q plot
    plt.figure()
    ax = plt.gca()
    stats.probplot(ncdata, dist='norm', fit=True, plot=ax, rvalue=True)
    plt.title(f'Quantile plot vs Normal Dist.     ({fullhash})')

    ########################
    # overplot the two cdfs
    my_dpi = 200
    w = 1250
    h = 1250
    plt.figure(figsize=(w/my_dpi, h/my_dpi), dpi=my_dpi)
    plt.plot(normdblevels, cumdist, normdblevels, ncdf) # compare two CDFs
    plt.title(f'CDF vs Normal Dist.     ({fullhash})')
    plt.show()




if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
