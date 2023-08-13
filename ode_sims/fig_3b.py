import PyDSTool
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
import matplotlib as mpl
from matplotlib.lines import Line2D

'''
This script runs the ode model for different membrane potentials and pH, and plots the
ATP synthase flux vs membrane potential (it generates Fig. 3B in the paper).
08/11/23
GCG
'''

params = {'axes.labelsize': 7,
           'axes.titlesize': 7,
          'legend.fontsize': 6,
           'xtick.labelsize': 7,
           'ytick.labelsize': 7,
            'figure.figsize':(1.8,1.8)}
mpl.rcParams.update(params)

Na = 6.02214e23
#Mito vol
vmito = 0.016e-15 #lt matrix
vims = 0.021e-15# IMS
vcube = 0.306e-15

#ADP and ATP concentration
cdm_i = 2 #10mM
ctm_i = 15 - cdm_i #mM
cdo_i = 0.1 #10mM
cto_i = 6.5 #mM

#ADP and ATP # of molecules
nr_dm_m = 0.45*0.8*6.02*cdm_i*vmito*1e20 #number of molecules Dm
nr_tm_m = 0.05*6.02*ctm_i*vmito*1e20
nr_do_i = 0.45*6.02*cdo_i*vims*1e20
nr_to_i =0.05*6.02*cto_i*vims*1e20
nr_to_c =0.05*6.02*cto_i*vcube*1e20

no1_ant = 16471 #number of ants
no_atp = 267 #number of atphase
k7 = 9.2 #0.74
k8 = 0.35 #0.05
k9 = 0.58 #0.37
k10 = 0.48#0.47
a12 = 100 #s-1
a65 = 3969
a21 = 40.0
a23 = 5.0 #muMs-1
a32 = 5e3
a16 = 146516
a61 = 33989
kp  = 1.0
n_porin = 6340

DSargs = PyDSTool.args(name='ex')

DSargs.pars = { 'k1_on':'4.0',#          #KDo
                'k1_off':'100', #        #KDo
                'k2_on': '6.4',#         #KTi
		        'k2_off':'4e4',#         #KTi
		        'k5_on':'0.4',           #KTo
	            'k5_off':'200.0',        #KTo
	         	'k6_on':'4.0',           #KDi
		        'k6_off':'4e4',          #KDi
		        'k7':k7,                 #kp
		        'k8':k8,                 #kcp
		        'k9':k9,                 #kt
		        'k10':k10,               #kd
		        'no_ant': no1_ant,
		        'a12':a12, #s-1
		        'a21':a21,#s-1
		        'a65':a65,             #s-1
		        'a56':'1e3',              #s-1,
		        'a16':a16,           #s-1
		        'a61':a61,            #s-1
		        'a25':'5.85e-30',         #s-1
		        'a52':'1e-20',            #s-1
		        'a54':'1e2',              #s-1
	         	'a45':'1e2',              #s-1
		        'a43':'2',                #uM-1s-1
		        'a34':'1e2',
		        'a32':a32,                #s-1
		        'a23':a23,                #uM-1s-1
		        'no1_atp': no_atp,        #number of atpases
		        'vims':vims,
		        'vcube':vcube,
		        'vmito':vmito,
		        'kp':kp,
		        'n_porin':n_porin,
		        'Na':Na,
		        'fa':'0.5',
		        'ac':nr_do_i,
                'am':nr_dm_m,
                'bm':nr_tm_m,
                'bo':nr_to_i ,

}

DSargs.varspecs = {'eo':'-a65*eo + a56*(no1_atp-eo-ei-h3eo-h3es-h3e)+ a16*ei - a61*eo',
                    'ei':'-a16*ei + a61*eo - a12*ei + a21*h3e',
                    'h3eo':'-a45*h3eo + a54*(no1_atp-eo-ei-h3eo-h3es-h3e) + a34*h3es - a43*h3eo*(1e6*am/(Na*vmito))',
                    'h3es':'a43*(1e6*am/(Na*vmito))*h3eo - a34*h3es + a23*h3e*(1e6*(bm)/(Na*vmito)) - a32*h3es',
                    'h3e':'-a23*h3e*(1e6*(bm)/(Na*vmito)) + a32*h3es - a25*h3e + a52*(no1_atp-eo-ei-h3eo-h3es-h3e) + a12*ei - a21*h3e'}

# initial conditions
DSargs.ics  = {'eo':no_atp,
                'ei':0,
                'h3eo':0,
                'h3es':0,
                'h3e':0}

DSargs.tdomain = [0,10]

#ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)    # an instance of the 'Generator' class.
ode  = PyDSTool.Generator.Radau_ODEsystem(DSargs)
traj = ode.compute('polarization')
pd   = traj.sample()


fig = plt.figure(1)
fig.subplots_adjust(right=0.98, left = 0.25, bottom =0.2, top = 0.95)

#color = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9']
#color = ['r','b','g','y']
color = ['0.2','0.4','0.6','0.8']

l_ph = [6.6, 7.0, 7.2,7.4]
for i,l in enumerate(l_ph):
    print(i,l, l - 7.6,color[i])
    a21 = 6.33e24*10**(-7.6*3)
    a65 = 1.58e25*10**(-l*3)

    for n,j in enumerate(np.arange(50,220,1)):
    #a21 = 6.33e24*10**(-7.6*3)
    #a65 = 1.58e25*10**(-j*3)
    #print(n,j, color[n])
        phim = j - 50 #mv
        a61 = 4.98e7*np.exp(-(3*phim)/(2*26.75))
        a16 = 1e2*np.exp((3*phim)/(2*26.75))
    #print(j-7.6, a65)


        ode.set( pars = {'a21':a21,
                         'a65':a65,
                         'a61':a61,
                         'a16':a16})

        tmp = ode.compute('pol%3i' % i).sample()    # or specify dt option to sample to sub-sample
        plt.plot(j,-(tmp['h3es'][-1]*a32*1e3)/(Na*vmito) + (tmp['h3e'][-1]*a23*ctm_i*0.05*1e6)/(Na*vmito) ,'.', c=color[i],markersize=1)

plt.xlim(40,220)
plt.ylim(-1,0.1)

plt.ylabel('ATP flux (mM/s)',labelpad =0.1)
plt.xlabel(r'$\Delta\psi$ (mV)',labelpad =0.1)

circ1 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.2",mec='0.2')
circ2 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.4",mec='0.4')
circ3 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.6",mec='0.6')
circ4 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.8",mec='0.8')
#circ5 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.5",mec='0.5')
#circ6 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.6",mec='0.6')
plt.legend((circ1, circ2,circ3,circ4), (r"$\Delta$pH -1",r"$\Delta$pH -0.6 ",r"$\Delta$pH - 0.4",r"$\Delta$pH - 0.2 "), numpoints=1,  frameon = False,bbox_to_anchor=(-0.08,0),loc="lower left")
#plt.tight_layout()

#plt.savefig('param_sweep_dph_phi7_6_flux_sgn_last.png', format = 'png',dpi=600,bbox_inches='tight')

plt.show()
