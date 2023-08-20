#
    ################################################
    #
    #    Now get to work...   #
    #
    #


    if SPACE == '2D':
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
            df.metadata.d['Computer Name'] = cto.PCNAME
            codeFolder = ''
            df.set_folders(pointsDataFolder,codeFolder) # also creates filename
            c1.fill(gt)  # randomize and compute traj's and costs
            df.metadata.d['grid info'] = f'{N}x{N} random grid, {N*N} pts.'
            gt.savePoints2D(df) #this will set the 'Research Question' to "RandomGridPointSet"
            print(f'Random pts saved to {df.name}')
            notes = f'generated random points file: {df.hashcode} {cto.N}x{cto.N}'
            logentry(df,notes)
            print(f'\n\n               your points data file hash is:',df.hashcode)

        elif OP_MODE == 'search': # search mode with rect grid or read-in points
            if gridtype == 'random': # we should read from file for repeatability
                # create the random points file reader
                dfr = bd.datafile('','','') #'' ok for reading
                dfr.set_folders('','') # '' ok for reading
                dfr.name = pointsFilename
                pointSourceHash = gt.readPoints2D(dfr)  #read in the set of random points
            c1.fill(gt) # calc trajectories and costs after points reading
            dfw = bd.datafile('2Dsearching','BH','simulation')
            dfw.metadata.d['Computer Name'] = cto.PCNAME
            dfw.set_folders(DataFolder,'')
            if gridtype=='random':
                dfw.metadata.d['Points Data Source'] = pointsHash
            else:
                dfw.metadata.d['Points Data Source'] = 'Regular Grid'
            q = input('Research Question for this search: ')
            dfw.metadata.d['Research Question'] = q
            print(f"RQ0: {dfw.metadata.d['Research Question']}")
            # instantiate a path:
            p = cto.path(gt,c1)
            path2, cmin = p.search(SEARCHT, dfile=dfw, nsamples=nsearch)
            print(f"RQ1: {dfw.metadata.d['Research Question']}")

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
            print(f"RQ2: {dfw.metadata.d['Research Question']}")
            #  keep a "log book"
            logentry(dfw,notes)
            # graph the optimal search result (best path)
            path2.plot(-1,notes)


    #    6D version
    #
    if SPACE == '6D':

        # create the datafile:
        df = bd.datafile('6Dsearching','BH','simulation')
        df.metadata.d['Computer Name'] = cto.PCNAME
        if gridtype == 'random':
            df.metadata.d['Random Grid'] = True
        else:
            df.metadata.d['Random Grid'] = False
        df.set_folders('/home/blake/Sync/Research/CalTrajOpt_RESULTS','')
        q = input('Research Question for this search: ')
        df.metadata.d['Research Question'] = q

        cto.pts = cto.setupPoints6D()   # just store points instead of cost matrix Cm
        #c1 = cto.Cm(df = df)  # cost matrix
        #c1.fill()

        # compute a path:
        #p = cto.path3d(gt,c1)
        p = cto.path3D()

        p.search(searchType,dfile=df,nsamples=nsearch,profiler=None)
        # search will close the datafile

        q = df.metadata.d['Research Question']
        notes = '{:}, grid: {:}, n={:}, {:}, cost: {:4.2f} ({:})'.format(searchType, cto.gridType, nsearch, SPACE, df.metadata.d['Min Cost'], cto.costtype)
        logentry(df,notes)
        ####

        print('Completed: see results at ',df.hashcode)

def logentry(df,notes):
    logdir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/'
    logfilename = 'work_logbook.txt'
    q = df.metadata.d['Research Question']
    if len(q)>0 and 'debug' not in q:  # skip junk files
        now = dt.datetime.now()
        dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
        logline = '{:}, {:}, {:},  RQ: {:}'.format(dtstring,df.hashcode, notes, q)
        with open(logdir+logfilename,'a') as f:
            print(logline,file=f)
            f.close()
            print('added log entry to: ',logdir+'work_logbook.txt')
    else:
        print(f'RQ: {q}')
        print(f'debugging detected. {df.hashcode} will not be logged to {logdir+logfilename}')
