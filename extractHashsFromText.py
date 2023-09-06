#

import re
import sys

if len(sys.argv) != 2:
    print ('usage: >python3 extractHashFromText  FILEname')
    quit()

f = open(sys.argv[1], 'r')

hashmatcher = r'(?<![0-9A-Fa-f])[0-9A-Fa-f]{8}(?![0-9A-Fa-f])' #thanks ChatGPT!
hmc = re.compile(hashmatcher)

hashset = set()

for line in f:
    results = hmc.findall(line)
    for result in results:
        hashset.add(str(result))
resultlist = sorted(list(hashset))
print(f'I found {len(hashset)} results:')

for r in resultlist:
    print(r)
