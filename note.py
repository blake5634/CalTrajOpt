import searchOpt as so
import brl_data.brl_data as bd

x = input('Entry for the worklog: ')
df = bd.datafile('','','')
df.metadata.d['Research Question'] = 'manual log entry'
df.hashcode = 'None    '
if len(x) > 0:
    so.logentry(df,x)
else:
    print(' ... log entry canceled')
