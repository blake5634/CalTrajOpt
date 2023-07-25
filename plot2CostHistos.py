#!/usr/bin/python3

#                                        Plot two cost histograms for comparison

import os
import c2to as cto
import sys
import numpy as np
import statistics as stats
from scipy.stats import norm
import brl_data.brl_data as bd
import matplotlib.pyplot as plt
import datetime as dt

def main(args):
    ## some basics
    #cto.configure()
    #gr = cto.grid2D(cto.N)
    #Cm = cto.Cm()
    #Cm.fill(gr)

    if len(sys.argv) != 3:
        bd.brl_error('usage: >python3 plot2CostHistos.py  hash1 hash2 ')
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

    ## select two files
    #print('Select the baseline (random) distribution to plot:')
    #hashstr1 = input('Enter first 4 characters of the hash:')
    #print('Select the 2nd (heuristic) distribution to plot:')
    #hashstr2 = input('Enter first 4 characters of the hash:')
    hashstr1 = sys.argv[1]
    hashstr2 = sys.argv[2]

    dfname1 = None
    dfname2 = None
    for f in csvfiles:
        if hashstr1 in f:
            dfname1 = f    # our desired file name
        if hashstr2 in f:
            dfname2 = f    # our desired file name
    if not dfname1:
        bd.brl_error('Somethings wrong with hash string '+hashstr1+' (not found)')
    if not dfname2:
        bd.brl_error('Somethings wrong with hash string '+hashstr2+' (not found)')

    #
    #  OK - now lets create a datafile and open it for reading
    #
    mds = [{},{}]
    hs = ['x','x']
    for ip, dfname in enumerate([dfname1, dfname2]):
        print('opening ', dfname)
        df = bd.datafile('', '','')  # open it with blank title info
        df.hashcode = dfname.split('_')[1]  # replace the newly generated hash code with the file of interest
        hs[ip]=df.hashcode
        df.set_folders(datadir,'')        # set these to wherever you want to open datafiles
        df.open('r',tname=datadir+'/'+dfname)
        df.metadata.polish()  # convert metadata from strings to useful types
        mds[ip] = df.metadata.d

    xs = [0,0]
    ys = [0,0]
    mus = [0,0]
    sds = [0,0]
    for ip, dfname in enumerate([dfname1, dfname2]):
        try:
            x =mds[ip]['CostHistogram_levels'][1:]
        except:
            cto.error('Please run getCostHisto in this data first')
        xs[ip] = x
        ys[ip] = mds[ip]['CostHistogram_values']
        mus[ip] = mds[ip]['CostMean']
        sds[ip] = mds[ip]['CostStDev']

    costtype = mds[0]['CostType']
    barwidth = 0.75*(max(x)-min(x))/len(x)
    colors = ['b','r']

    #
    # do the plotting
    #

    my_dpi = 200
    w = 1800
    h = 650
    plt.figure(figsize=(w/my_dpi, h/my_dpi), dpi=my_dpi)
    maxx = 0
    for ip, dfname in enumerate([dfname1, dfname2]):
        if max(xs[ip]) > maxx:
            maxx = max(xs[ip])
        plt.bar(xs[ip],ys[ip],width=barwidth,color=colors[ip])
        #overplot normal distribution
        curve = norm.pdf(xs[ip],mus[ip],sds[ip])
        scale = max(ys[ip])/max(curve)
        for i in range(len(curve)):
            curve[i] *= scale
            #c2[i] *= scale
        plt.plot(xs[ip],curve)
    #round up maxx by 1000
    maxx = (maxx//1000 + 1)*1000
    # plot relative to 0 for visual comparison
    plt.xlim([0,maxx])
    plt.title('{:} cost distributions'.format(costtype))
    plt.xlabel('Cost      ({:},{:})'.format(hs[0],hs[1]))
    plt.ylabel('# paths')
    plt.tight_layout() # make sure stuff shows with custom dimensions
    figure = plt.gcf()
    plt.show()

    template = '______________'+hs[0]+'_'+hs[1]+'.png'
    print('filename template: ',template)
    nroot = input('enter name root: (<enter> to not save) ')
    if len(nroot)>0:
        h0 = hs[0]
        h1 = hs[1]
        nroot += '_' # separate the hash
        imgdir = datadir+'/writing/'
        imgname = nroot + hs[0]+'_'+hs[1]+'.png'
        idpi = 200
        cto.plotSave(figure, idpi, imgdir, imgname)

    else:
        print('plot image NOT saved')



if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
