#!/usr/bin/python3

import glob
import os
import c2to as cto
import sys
import brl_data.brl_data as bd

def main(args):
    #####################################################################
    #
    # get all the files
    #
    datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    # thanks stackoverflow
    files = list(glob.glob(datadir + '/' + "*"))
    files.sort(key=lambda x: os.path.getmtime(x),reverse=True) # newest first
    filenameroots = []
    if len(files) ==0:
        cto.error('No brl_data files found')
    for f in files:
        #print('found: ',f)
        if '.csv' in f :
            filenameroots.append(str(f).replace('.csv',''))
        if  '_meta.json' in f:
            filenameroots.append(str(f).replace('_meta.json',''))
    # dict.fromkeys better than set because preserves order (thanks google)
    filenameroots = list(dict.fromkeys(filenameroots)) # elim dupes


    #
    #  OK - now lets check all metadatas
    #
    remlist = []
    remlistCSV = []
    for fnr in filenameroots:
        mdname = fnr+'_meta.json'
        if(os.path.exists(mdname)):  # if we can find a metadata file
            df = bd.datafile('', '','')  # open it with blank title info
            df.set_folders(datadir,'')        # set these to wherever you want to open datafiles
            df.open(mode='r',tname=mdname)
            #df.metadata.polish()  # convert metadata from strings to useful types
            df.close()
            # RQ key came in two forms: with and without a space
            try:
                q = df.metadata.d['Research Question']
                if len(q)==0 or 'debug' in q.lower(): # there is no real RQ
                    remlist.append(fnr+'_meta.json')
                    remlist.append(fnr+'.csv')
            except:
                pass
            try:
                q = df.metadata.d['ResearchQuestion']
                if len(q)==0 or 'debug' in q.lower(): # there is no real RQ
                    remlist.append(fnr+'_meta.json')
                    remlist.append(fnr+'.csv')
            except:
                pass
        else:
            print('non existent metadata file: ',fnr)
            if 'ties' not in fnr:  # we have special csv files for ties data
                remlistCSV.append(fnr+'.csv')

    if len(remlist) > 0:
        print('\n\nPlanning to remove the following file names (.jason and .csv):')
        for fn in remlist:
            print(fn)

        x = input('\n\n          OK to delete these files?? ... (y/N)')
        if x.lower() == 'y':
            for fn in remlist:
                print(' ... removing ',fn)
                os.remove(fn)
        else:
            print('removing canceled')

    if len(remlistCSV) > 0:
        print('\n\nPlanning to remove the following .csv files which have no metadata:')
        for fn in remlistCSV:
            print(fn)

        x = input('\n\n          OK to delete these files?? ... (y/N)')

        if x.lower() == 'y':
            for fn in remlistCSV:
                print(' ... removing ',fn)
                os.remove(fn)
        else:
            print('removing canceled')

    print('.. Done')






if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
