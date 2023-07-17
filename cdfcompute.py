import scipy.stats as stat
import numpy as np

# 9.99 dummy value for future
x = [3.86, 3.43,3.85, 3.34,3.86,2.9,4.84,4.36,6.5,5.44,4.36, 4.29]
y = []
for v in x:
   y.append(100.0*(1.0 - stat.norm.cdf(v)))
for i,v in enumerate(x):
   print('{:6.3f}sig = {:14.8f} pct'.format(v,y[i]))
