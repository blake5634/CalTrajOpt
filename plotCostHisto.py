#!/usr/bin/python3

#                                        Plot cost histogram

import os
import c2to as cto
import sys
import statistics as stats
from scipy.stats import norm
import brl_data.brl_data as bd
import matplotlib.pyplot as plt

def main(args):
    # some basics
    cto.configure()
    gr = cto.grid2D(cto.N)
    Cm = cto.Cm()
    Cm.fill(gr)
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

    x = df.metadata.d['CostHistogram_levels'][1:]
    y = df.metadata.d['CostHistogram_values']
    npts = sum(y)
    mu,sd = 17.86, 2.34
    plt.bar(x,y,width=0.5,color='b')
    curve = norm.pdf(x,mu,sd)
    scale = npts
    for i in range(len(curve)):
        curve[i] *= scale
    plt.plot(x,curve)
    plt.show()


    #p.plot(idx)





if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
