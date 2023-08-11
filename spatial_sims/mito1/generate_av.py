import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter

'''
This script generates the average over 10 initial conditions
GCG
8.11.23
'''

seed_number = 10 #number of seeds
nv = 30
npoints = int(300000)#int(1e5) #time points available
arr = np.zeros((nv,npoints,seed_number))
av_ar = np.zeros((nv,npoints)) #average traces
ar_std = np.zeros((nv,npoints))
time = np.zeros((nv, npoints)) #time

li = ['/D.outer_membrane_final.dat','/D.inner_membrane_final.dat','/T.inner_membrane_final.dat','/T.outer_membrane_final.dat','/L.inner_membrane_final.dat',
'/DL.World.dat','/LD.World.dat','/TL.World.dat','/LT.World.dat','/DLT.World.dat',#5 first
'/TLD.World.dat','/DLD.World.dat','/DLDp.World.dat','/TLT.World.dat','/TLTp.World.dat',#10 first
'/Eo.World.dat','/Ei.World.dat','/H3Eo.World.dat', '/H3E.World.dat','/EH3.World.dat',#15
'/H3ES.World.dat','/atp_prod.World.dat','/atp_dis.World.dat','/prod.World.dat','/counter_prod.World.dat',#20
'/unprod_d.World.dat','/unprod_dp.World.dat','/exp_t.World.dat','/imp_t.World.dat','/T.Cube.dat']#25


for i in range(seed_number):
    print(i)
    if i>8:
        for s,j in enumerate(li):
            var = np.genfromtxt("./react_data/seed_000"+str(i+1)+j, dtype = float)
            arr[s,:,i]= var[:npoints,1]
            time[s,:] = var[:npoints,0]

    else:
        for s,j in enumerate(li):
            var = np.genfromtxt("./react_data/seed_0000"+str(i+1)+j, dtype = float)
            arr[s,:,i]= var[:npoints,1]
            time[s,:] = var[:npoints,0]

arr[29,:,:] = arr[29,:,:]-arr[3,:,:] #T cube
arr[0,:,:] = arr[0,:,:]-arr[1,:,:]      #D IMS
arr[3,:,:] = arr[3,:,:]-arr[2,:,:]     #T IMS
arr[21,:,:] = arr[21,:,:]-arr[22,:,:]
arr[23,:,:] = arr[23,:,:]-arr[24,:,:]
arr[27,:,:] = arr[27,:,:]-arr[28,:,:]

av_ar = np.average(arr, axis=2)
ar_std = np.std(arr, axis=2)

the_filename = 'av_10c_1e8_m1'
with open(the_filename, 'wb') as f:#
	pickle.dump(av_ar, f)#

#np.savetxt('time_10c_1e8',time)
np.savetxt('t_cyto_10c_1e8',av_ar[29,:])
np.savetxt('d_ims_10c_1e8',av_ar[0,:])
np.savetxt('t_ims_10c_1e8',av_ar[3,:])
np.savetxt('d_matrix_10c_1e8',av_ar[1,:])
np.savetxt('t_matrix_10c_1e8',av_ar[2,:])
np.savetxt('time_10c_1e8',time[0,:])
