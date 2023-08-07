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
    ## some basics
    #cto.configure()
    #gr = cto.grid2D(cto.N)
    #Cm = cto.Cm()
    #Cm.fill(gr)
    if len(sys.argv) != 2:
        cto.error('usage: >p3 plotCostHisto <hashcode>')
    hashstr = sys.argv[1]

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
    df.hashcode = dfname.split('_')[1]  # replace the newly generated hash code with the file of interest
    print('Our hashcode is: ',df.hashcode)
    df.set_folders(datadir,'')        # set these to wherever you want to open datafiles
    df.open('r',tname=datadir+'/'+dfname)
    df.metadata.polish()  # convert metadata from strings to useful types

    try:
        x = df.metadata.d['CostHistogram_levels'][1:]
    except:
        cto.error('Please run getCostHisto in this data first')
    y = df.metadata.d['CostHistogram_values']
    npts = sum(y)
    barwidth = 0.75*(max(x)-min(x))/len(x)
    mu = df.metadata.d['CostMean']
    sd = df.metadata.d['CostStDev']
    ct = df.metadata.d['CostType']
    hc = df.hashcode

    my_dpi = 200
    w = 1250
    h = 1250
    plt.figure(figsize=(w/my_dpi, h/my_dpi), dpi=my_dpi)
    plt.bar(x,y,width=barwidth,color='b')
    plt.title('{:} cost distribution, n={:} samples'.format(ct,npts))
    plt.xlabel('Cost      ({:})'.format(hc))
    plt.ylabel('# paths')
    plt.tight_layout() # make sure stuff shows with custom dimensions
    curve = norm.pdf(x,mu,sd)
    #c2 = norm.pdf(x,18.2,sd)
    txtx = mu+1.3*sd
    txty = max(y)*0.8
    plt.text(txtx,txty,'mu: {:4.2f} sd: {:4.2f}'.format(mu,sd))
    # chatGPT: max_value / (1 / (σ * sqrt(2π))).
    scale = max(y)/max(curve)
    for i in range(len(curve)):
        curve[i] *= scale
        #c2[i] *= scale
    plt.plot(x,curve)
    fig = plt.gcf() #save the figure in RAM
    plt.show()

    template = '______________'+df.hashcode+'.png'
    print('filename template: ',template)
    nroot = input('enter name root: (<enter> to not save) ')
    if len(nroot)>0:
        nroot += '_' # separate the hash
        imgdir = datadir+'/writing/'
        imgname = nroot + df.hashcode
        cto.plotSave(fig, my_dpi, imgdir, imgname)

    else:
        print('plot image NOT saved')



    #p.plot(idx)





if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
