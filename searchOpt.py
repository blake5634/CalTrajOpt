#!/usr/bin/python3

#import numpy as np
import c2to as cto
import sys
import brl_data.brl_data as bd
import datetime as dt


#  Usage:
#
# 1) Generate points data file (for random grid re-use)
#
#    >p3 searchOpt.py generate
#
# 2) Search a stored random grid
#
#    >p3 searchOpt.py hhhhhh  (hash code of grid points file)
#
# 3) Search a rectangular grid
#
#    >p3 searchOpt.py
#

def main(args):
    cto.configure()
    ###################################################
    #    Process command line
    #
    # this code can do two things: generate grid, search grid.
    #  use 'generate' CL arg to generate a fresh random grid
    #  use hashcode to specify random points file
    if len(args) == 2:   # >p3 searchOpt.py  [generate | hashcode ]
        gridtype = 'random'   # 2 arg version always 'random'
        if args[1].startswith('gen'):
            OP_MODE = 'generate'
        else: # we got a hashcode
            OP_MODE = 'search'
            pointsHash = args[1]
            print(f'Reading random points from {pointsHash}')
    if len(args) == 1:   # >p3 searchOpt.py  (search rect grid)
        OP_MODE = 'search'
        gridtype = 'rectangular'  # 1 arg version never random grid

    ###################################################  checks
    if OP_MODE not in ['generate', 'search']:
        cto.error('searchOpt.py: illegal POINTS MODE')

    if OP_MODE == 'generate' and gridtype != 'random':
        cto.error('SearchOpt: Cannot generate points file unless grid is random')
    ###################################################

    DataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    PointsFolder = '/PointSetsRandom'
    pointsDataFolder = DataFolder + PointsFolder
    #pointsDataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom'
    #pointsFilename = pointsDataFolder + '/'+pointsDataName  # full path of points file

    if OP_MODE == 'search' and gridtype == 'random':  # we need to read points file
        myf = bd.finder()
        myf.set_dirs([DataFolder+PointsFolder])
        keys = [pointsHash, '.csv']            # look for requested rand points file
        fns = myf.findh(keys)
        pointsDataFolder = fns[0][0]  # should be same as .set_dirs above
        pointsDataName   = fns[0][1]
        pointsFilename = pointsDataFolder + '/'+pointsDataName  # full path of points file

    ##########################################################################
    #
    #    Configure the job
    #
    cto.N = 3
    cto.M = cto.N*cto.N
    N = cto.N

    #
    #   Choose search type
    #
    #SEARCHT = 'heuristic search' # greedy nearest neighbor
    #SEARCHT = 'exhaustive'   # enumerate all paths (formerly 'brute force')
    #SEARCHT = 'sampling search' # nsearch random paths
    SEARCHT = 'multi heuristic' # repeated heuristic search all starting pts
    #
    #   Choose search size
    #
    #nsearch = int(np.math.factorial(9)* 0.10)  # 10% of 3x3
    #nsearch = 1000000  # 1M
    nsearch = 4*N*N   # 4 searches from each starting pt
    print(f'Test: {nsearch} searches planned')
    #nsearch = 10

    #
    #   Choose cost type
    #
    cto.costtype = 'time'
    cto.costtype = 'energy'
    cto.NPC = 30   #  # of simulation points in 0-dt time interval
    #
    ##########################################################################

    #
    #    Now get to work...   #
    #
    # create grid
    gt = cto.grid2D(cto.N)
    # create cost matrix
    c1 = cto.Cm()
    if gridtype == 'random':
        gt.randgrid = True
        c1.set_GridRandomize()  # select random instead of grid

    if OP_MODE == 'generate':
        if gridtype != 'random':
            cto.error('Somethings wrong with search config: generate not random')
        df = bd.datafile(f'2DrandomGrid{N}x{N}PointSet','BH','simulation')
        codeFolder = ''
        df.set_folders(pointsDataFolder,codeFolder) # also creates filename
        c1.fill(gt)  # randomize
        df.metadata.d['grid info'] = f'{N}x{N} random grid, {N*N} pts.'
        gt.savePoints2D(df) #this will set the 'Research Question' to "RandomGridPointSet"
        print(f'Random pts saved to {df.name}')
        notes = f'generated random points file: {df.hashcode} {cto.N}x{cto.N}'
        logentry(df,notes)
        print(f'\n\n               your points data file hash is:',df.hashcode)

    elif OP_MODE == 'search': # search mode with grid or read-in points
        if gridtype == 'random': # we should read from file for repeatability
            # create the random points file reader
            dfr = bd.datafile('','','') #'' ok for reading
            dfr.set_folders('','') # '' ok for reading
            dfr.name = pointsFilename
            pointSourceHash = gt.readPoints2D(dfr)  #read in the set of random points
        c1.fill(gt) # calc trajectories and costs after points reading
        dfw = bd.datafile('2Dsearching','BH','simulation')
        dfw.set_folders(DataFolder,'')
        if gridtype=='random':
            dfw.metadata.d['Points Data Source'] = pointsHash
        else:
            dfw.metadata.d['Points Data Source'] = 'Regular Grid'
        q = input('Research Question for this search:')
        dfw.metadata.d['Research Question'] = q
        # instantiate a path:
        p = cto.path(gt,c1)
        path2, cmin = p.search(SEARCHT, dfile=dfw, nsamples=nsearch)
        print('Optimal path returned: (tra)', path2.path)
        print('Optimal path returned: (idx)', path2.idxpath)
        # is it a valid path?
        #p.check()
        if gridtype=='random':
            notes = f"Search Result: Rand-grid({pointSourceHash}), {SEARCHT}, cost: {cmin:.1f} ({cto.costtype})"
            print(f'\n\n               your search results file hash is: {dfw.hashcode} using grid {pointSourceHash}.')
        else:
            notes = f"Search Result: {gridtype} grid, {SEARCHT}, cost: {cmin:.1f} ({cto.costtype})"
            print(f'\n\n               your search results file hash is: {dfw.hashcode}.')
        #  keep a "log book"
        logentry(dfw,notes)
        # graph the optimal search result (best path)
        path2.plot(-1,notes)


def logentry(df,notes):
    logdir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/'
    logfilename = 'work_logbook.txt'
    q = df.metadata.d['Research Question']
    if len(q)>0 and 'debug' not in q:  # skip junk files
        now = dt.datetime.now()
        dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
        logentry = '{:}, {:}, {:},  RQ: {:}'.format(dtstring,df.hashcode, notes, q)
        f =open(logdir+logfilename,'a')
        print(logentry,file=f)
        f.close()
    else:
        print(f'debugging detected. {df.hashcode} will not be logged to {logdir+logfilename}')

if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
