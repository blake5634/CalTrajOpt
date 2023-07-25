
#  convert between 6D int coordiates and point index

N=3

def getidx(v):
    return v[0]*N**5+v[1]*N**4+v[2]*N**3+v[3]*N**2+v[4]*N+v[5]

def getcoord(idx):
    v = [0,0,0,0,0,0]
    r = idx
    for i in range(6):
        v[5-i] = r%N
        r = (r-v[5-i])//N
    return v

for i in [0,700,701,702,703,727,728]:
    assert (getidx(getcoord(i))== i), 'assertion '+str(i)+'failed!'
    print('{:3} '.format(i),getcoord(i))
