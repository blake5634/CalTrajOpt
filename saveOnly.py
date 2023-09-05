#!/usr/bin/python3

import glob
import os
import c2to as cto
import sys
import brl_data.brl_data as bd
import datetime as dt

def main(args):

    #  all the dirs we might find data files
    #dirs = ['/home/blake/Sync/Research/CalTrajOpt_RESULTS',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D-Round2',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data/Gold',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing'
            #]

    # for testing
    dirs = ['/home/blake/Ptmp/CalTrajOpt/testing/folder1',
            '/home/blake/Ptmp/CalTrajOpt/testing/folder2',
            '/home/blake/Ptmp/CalTrajOpt/testing/folder3'
            ]

    flags = ['nothing']

    if len(args) != 2:
        bd.brl_error('usage: >saveOnly.py   URIfile')


    # 0) collect the URLs to save
    uriFile = args[1]
    try:
        f = open(uriFile,'r')
    except:
        bd.brl_error(f'URI file: {uriFile} not found')
    URIs2Save = []
    for l in f:
        l = l.strip()
        if len(l) < 6 or len(l)>8:
            bd.brl_error(f'URI:{l} is less than minimum 6 characters or > 8')
        for c in l:
            if c not in 'abcdef0123456789':
                bd.brl_error(f'illegal URI charcter: {c} in {l}')
        URIs2Save.append(l)

    # 1) find all files NOT in the URIs2Save list
    fdList = []
    hset = set()
    for d in dirs:
        fs = os.listdir(d) # file names
        for f in fs:
            include = True
            thishash = bd.getHashFromFilename(f)
            for u in URIs2Save:
                print(f'file {f}: comparing {u} vs {thishash}')
                if u == thishash:
                    include = False # were saving these!
            print(f'           include: ',include)
            if include:
                print('             appending: ',thishash)
                fdList.append([d,f])
                hset.add(thishash)

    hlist = list(hset)
    # 2) purge those files

    hrem = purgeFilesbyHash(hlist,dirs)
    for h in hrem:
        log_record_deletions(h,'file')


    # 3) find all log entries NOT having URIs2Save in them

    # 4) purge the log entries
    hrem = purgeLogsbyHash(hlist,flags=flags)
    for h in hrem:
        log_record_deletions(h,'log entry')
    quit()


def purgeFilesbyHash(hlist,dirs,flags=['None']):
    hashesRemoved = set()
    myf = bd.finder()
    myf.set_dirs(dirs)
    fds = [] # 'fd' is a pair (dir,fname)
    for h in hlist:
        keys = [h]
        fds += myf.findh(keys) # get all files with this hash
    print(f' I found {len(fds)} files matching {hlist}:')
    i=0
    for fd in fds:
        i+=1
        print(f'    {i}: [{fd[0]}/{fd[1]}]')
    if 'listOnly' in flags:
        return []
    x = input('\n\n          OK to delete these files?? ... (y/N)')
    if x.lower() == 'y':
        for fd in fds:
            fn = fd[1]
            print(' ... removing ',fn)
            hashesRemoved.add(bd.getHashFromFilename(fn))
            os.remove(fd[0]+'/'+fd[1])
        return list(hashesRemoved)
    else:
        print('removing canceled')
        return []


def purgeLogsbyHash(hlist,flags=['None']):
    print(f'\nChecking  {hlist} in log files.')
    # now purge entries from the logs
    wlog = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/work_logbook.txt'
    ilog = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/image_log.txt'
    logs = [wlog, ilog]
    for lf in logs:
        justname = lf.split('/')[-1]
        f = open(lf,'r')
        keeplines = []
        deletelines = []
        for line in f:
            save = True
            for h in hlist:
                if h in line:
                    save = False
            if save:
                keeplines.append(line)
            else:
                deletelines.append(line)
        if 'listOnly' in flags:
            pass
            #print(f'\n    Matched {len(deletelines)} {hlist} entries from {justname}')
        else:
            print(f'\n    Deleting {len(deletelines)} {hlist} entries from {justname}')
        f.close()

        # explain what you will do
        print(f'{justname}:')
        i=0
        for l in deletelines:
            i+=1
            l = l.strip()
            print(f'    {i}: [{l}]')

        # confirm y/N ??
        if 'listOnly' not in flags:
            x = input(f'\n   OK to delete {len(deletelines)} entries from {justname}? (y/N): ')
            if x.lower() == 'y':
                # now modify the file
                f = open(lf,'w')
                for l in keeplines:
                    print(l.strip(),file=f)
                f.close()
                print(f' {len(deletelines)} lines purged from {justname}.')
            else:
                print(' purging logs canceled.')
    if 'listOnly' in flags:
        print('\n\n             List only, no actions were taken.')
    return [] # we don't log hashes deleted from log files

def log_record_deletions(hash, dtype):
    logfile = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/deletions.txt'
    f = open(logfile, 'a')
    now = dt.datetime.now()
    dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
    logline = f'{dtstring}: deleted {dtype} for {hash}'
    print(logline,file=f)
    f.close()

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
        print(f'debugging detected. {df.hashcode} will not be logged to {logdir+logfilename}')



if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
