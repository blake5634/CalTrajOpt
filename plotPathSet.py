#!/usr/bin/python3

import os
import c2to as cto
import sys
import brl_data.brl_data as bd

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
    print(df.metadata)

    print('')
    print('first rows of data are:')
    df.print_header()
    n=0

    #
    #   selection of paths:
    #
    c1 = 22
    c2 = 999  # 1st half of 1st quartile

    selectedpaths = []
    for row in df.reader:  # this is set up by the line df.open('r',tname=dfname)
        df.print_row(row)  # automatically formats the row according to types in metadata
        c = row[9]
        if type(c) != type(3.14):
            cto.error('cost type error')
        if c1 < c and c < c2:
            newp = cto.path(gr,Cm)   # add this path b/c its in the cost range
            pidx = row[:-1]  # path as node indices (drop cost)
            for k,pi in enumerate(pidx):  # make a path from this record
                if k==0:
                    next
                i1,j1 = cto.idx2ij(pidx[k-1])
                i2,j2 = cto.idx2ij(pidx[k])
                p1 = cto.point2D(i1,j1)
                p2 = cto.point2D(i2,j2)
                #print(p1,' --> ',p2)
                tr = cto.trajectory2D(p1,p2)
                tr.constrain_A() # make it right speed for AMAX
                newp.path.append(tr)
            if len(newp.path) != 9:
                cto.error('Path is wrong length: {:}/9'.format(len(newp.path)))
            selectedpaths.append(newp)

    # graph the path
    notes = 'Slowest paths overplotted'
    if len(selectedpaths) < 1:
        cto.error('select found no paths')
    else:
        print('selected {:} paths'.format(len(selectedpaths)))
    fig = selectedpaths[0].plotSetup(notes)  # setup
    for p in selectedpaths:
        p.plotOnePath(fig)
    p.plotDone()

    #p.plot(idx)





if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
