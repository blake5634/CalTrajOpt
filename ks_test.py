import math
import numpy as np
import scipy.stats as stats
import random
import sys
import os
import brl_data.brl_data as bd

#   Apply the Kolmogorov-Smirnov test of normality to cost histogram

testMode = False  #mess up the data to test the stat

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

    tvalues = []

    for row in df.reader:
        tvalues.append(float(row[-1]))

    mu = np.mean(tvalues)
    sd = np.std(tvalues)
    n = len(tvalues)

    print('orig mean, sd: ',mu,sd, 'n=',n)

    costvalues = []
    for d in tvalues:
        costvalues.append((d-mu)/sd)
        #print(' c =',row[-1], type(row[-1]))

    print('norm mean,sd: ',np.mean(costvalues),np.std(costvalues))

    if testMode:
        addptscount = int(len(costvalues)*1.5)
        for i in range(addptscount):
            costvalues.append(random.uniform(1000,2000))


    ##     perform Kolmogorov-Smirnov test for normality
    #k = kstest(random.uniform(1000,2000), norm.cdf)
    k = stats.ks_1samp(costvalues,stats.norm.cdf)
    ChiSq = k.statistic
    pval = k.pvalue


    print('Koloogorov-Smirnov results: ChiSq: {:8.5f}  pval = {:8.5f}, n={:}'.format(ChiSq, pval, len(costvalues)))

if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
