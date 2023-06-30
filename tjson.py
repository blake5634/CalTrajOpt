#!/usr/bin/python3

#
#  little very simple tests of python json serialization
#
import json

testfl = 3.14159
testint = 7
teststr = 'hello'
testlist = ['once', 'upon', 'a', 'time']
testlist2 = [5.1,6.2,7.3,8.4]

class t:
    def __init__(self):
        self.d = {}
        self.d['testfloat']  = testfl
        self.d['testint']    = testint
        self.d['teststring'] = teststr
        self.d['testlist1']  = testlist
        self.d['testlist2']  = testlist2

testobj = t()

f = open('jtestdata.json','w')
json.dump(testobj.d,f,indent=3)
f.close()


f = open('jtestdata.json','r')
d2 = json.load(f)
print('{:20}  | {:35} | {:}'.format('Key','Value','Type'))
print('-----------------------------------------------------------------------------------')
for k in d2.keys():
    v = d2[k]
    print('{:20}  | {:35} | {:}'.format(k,str(v),type(v)))
f.close()


