import pandas as pd
import  matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

'''
This script generates Fig 8C from the paper.
8.11.23
GCG
'''

params = {'axes.labelsize': 8,
           'axes.titlesize': 8,
          'legend.fontsize': 8,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
            'figure.figsize': (2,2)}
mpl.rcParams.update(params)

data = pd.read_csv('rate_all.csv')


lr_m1 = np.genfromtxt("./linear_regression_all_M1.txt")
lr_m12 = np.genfromtxt("./linear_regression_all_M12.txt")
lr_m19 = np.genfromtxt("./linear_regression_all_M19.txt")
lr_m7 = np.genfromtxt("./linear_regression_all_M7.txt")
lr_m14 = np.genfromtxt("./linear_regression_all_M14.txt")
lr_m4 = np.genfromtxt("./linear_regression_all_M4.txt")
lr_m6 = np.genfromtxt("./linear_regression_all_M6.txt")
lr_m10 = np.genfromtxt("./linear_regression_all_M10.txt")
lr_m11 = np.genfromtxt("./linear_regression_all_M11.txt")

print(np.mean(lr_m6[:,0]),np.mean(lr_m10[:,0]),np.mean(lr_m11[:,0]))

l1 = []
l1.extend(lr_m1[:,0])
l1.extend(lr_m6[:,0])
l1.extend(lr_m11[:,0])
l1.extend(lr_m12[:,0])
l1.extend(lr_m19[:,0])

l2 =[]
l2.extend(lr_m4[:,0])
l2.extend(lr_m7[:,0])
l2.extend(lr_m10[:,0])
l2.extend(lr_m14[:,0])


l = [np.mean(lr_m1[:,0]),np.mean(lr_m12[:,0]),np.mean(lr_m19[:,0]),np.mean(lr_m7[:,0]),np.mean(lr_m14[:,0]),np.mean(lr_m4[:,0]),np.mean(lr_m6[:,0]),np.mean(lr_m11[:,0]),np.mean(lr_m10[:,0])]
sd = [np.std(lr_m1[:,0]),np.std(lr_m12[:,0]),np.std(lr_m19[:,0]),np.std(lr_m7[:,0]),np.std(lr_m14[:,0]),np.std(lr_m4[:,0]),np.std(lr_m6[:,0]),np.std(lr_m11[:,0]),np.std(lr_m10[:,0])]

fig, ax = plt.subplots(figsize=(2, 2))
fig.subplots_adjust(right=0.9, left = 0.25, bottom =0.15, top = 0.85)

plt.errorbar(data['cluster'],data['Slope_av'] ,yerr= sd,fmt ='o', color='r',ms=5)

#plt.xlabel(r'Cluster')
plt.xlim(-0.1,1.1)
plt.ylim(0,40)
plt.ylabel('ATPs/msec in the Cytosol')
plt.xticks([0,1],['Globular','Elongated'])
#ax.set_xticklabels(,fontsize=9)
#plt.savefig('glob_atp_rate.png',dpi=600)

plt.show()
