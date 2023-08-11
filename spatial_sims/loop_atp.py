from scipy import stats
from matplotlib.lines import Line2D
import numpy as np

'''
This script was used to calculate the slope of the curve # ATP molecules vs time, and gives the rate of ATP production.
GCG
8.11.23
'''

seed_number = 10
f = open('linear_regression_all_M1.txt','w')

for i in range(seed_number):
    print(i+1)
    if i>8:
        tc = np.genfromtxt("./react_data/seed_000"+str(i+1)+"/T.Cube.dat", dtype = float)
        to = np.genfromtxt("./react_data/seed_000"+str(i+1)+"/T.outer_membrane_final.dat", dtype = float)
    else:
        tc = np.genfromtxt("./react_data/seed_0000"+str(i+1)+"/T.Cube.dat", dtype = float)
        to = np.genfromtxt("./react_data/seed_0000"+str(i+1)+"/T.outer_membrane_final.dat", dtype = float)

    slope, intercept, r_value, p_value, std_err= stats.linregress(tc[25000:100000,0]*1e3,tc[25000:100000,1]-to[25000:100000,1])
    print(p_value,tc[25000:100000,0]*1e3)
    f.write(str(slope))
    f.write('\t')
    f.write(str(r_value))
    f.write('\t')
    f.write(str(std_err))
    f.write('\n')
