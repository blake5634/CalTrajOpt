#!/usr/bin/python3

#                                        Plot cost histogram

import os
import c2to as cto
import sys
import numpy as np
import statistics as stats
from scipy.stats import norm
import brl_data.brl_data as bd
import matplotlib.pyplot as plt

def main(args):

    ##
    #
    #      create synthetic gaussian data file and then see if it hisograms correctly
    #           with plotCosthistogram.py
    #
    ##
    #
    #  OK - now lets create a datafile and open it for writing
    #
    npts = 100000
    N = 3
    expo = 2  #  2 or 6
    mu  = 30
    sig = 3


    datadir = '/home/blake/Ptmp/CalTrajOpt'


    print('opening new df')
    df = bd.datafile('normtest', 'BH','simulation')  # open it with blank title info
    df.set_folders(datadir,'')        # set these to wherever you want to open datafiles

    itype = str(type(5))
    ftype = str(type(3.1415))
    #tps = [itype]*(N**expo)      # path point seq
    tps = []
    for i in range(N**expo):
        tps.append(itype)
    tps.append(ftype) # the path cost's type

    names = []
    for i in range(N**expo):
        names.append('{:}'.format(i))
    names.append('Cost')
    df.metadata.d['Ncols'] = len(names)
    df.metadata.d['Types'] = tps
    df.metadata.d['Names'] = names
    df.metadata.d['CostType'] = 'time'
    df.metadata.d['SearchType'] = 'sampling'
    df.metadata.d['#samples'] = npts
    df.metadata.d['Research Question'] = ''  # triggers delete with cleanup script

    df.open('w')
    dfname = df.name
    print('new name: ',df.name)
    df.hashcode = dfname.split('_')[1]  # replace the newly generated hash code with the file of interest
    thash = df.hashcode # save it
    print('Our hashcode is: ',df.hashcode)

    costs = []
    for r in range(npts):
        trow = []
        for i in range(N**expo):
            trow.append(1)   # place holder data
        cost = np.random.normal(mu,sig)
        costs.append(cost)
        trow.append(cost)  # cost is normal random variable
        df.write(trow)

    df.metadata.d['Min Cost'] = min(costs)
    df.metadata.d['Max Cost'] = max(costs)
    df.close()

    print('File saved: hash: ',df.hashcode)
    quit()

    # after this runs
    #  1) getCostHisto
    #  2) plotCostHisto




if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
