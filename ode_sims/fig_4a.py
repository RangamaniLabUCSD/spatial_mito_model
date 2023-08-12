import PyDSTool
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
import matplotlib as mpl
from matplotlib.lines import Line2D

'''
ANT flux at steady state. This script generates the figure 4a in the paper
08.11.23
GCG
'''
params = {'axes.labelsize': 7,
           'axes.titlesize': 7,
          'legend.fontsize': 7,
           'xtick.labelsize': 7,
           'ytick.labelsize': 7,
            'figure.figsize': (1.8,1.8)}
mpl.rcParams.update(params)

Na = 6.02214e23

vmito = 0.016e-15 # matrix vol
vims =  0.021e-15 # IMS vol
vcube = 0.306e-15 # cube vol

#ADP & ATP concentrations
cdm_i = 2 #10mM
ctm_i = 13 #mM
cdo_i = 0.1 #10mM
cto_i = 6.5 #mM

#ADP & ATP number of molecules
nr_dm_i = 0.45*0.8*6.02*cdm_i*vmito*1e20 #number of molecules Dm
nr_tm_i = 0.05*6.02*ctm_i*vmito*1e20
#print(dm_i,tm_i)
nr_do_i = 0.45*6.02*cdo_i*vims*1e20
nr_to_i =0.05*6.02*ctm_i*vims*1e20

no1_ant = 16471 #number of atphase
no_atp = 267 #number of atphase
k7 = 9.20 #0.74
k8 = 0.35 #0.05
k9 = 0.58 #0.37
k10 = 0.48#0.47
a21 = 40.0
a23 = 5.0
a32 = 5e3
a12 = 25 #s-1
a65 = 3969
kp  = 0.5
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
		        'vims':vims,
		        'vcube':vcube,
		        'vmito':vmito,
		        'Na':Na,
		        'fa':'0.5',
		        'ac':nr_do_i,
                'am':nr_dm_i,
                'bc':nr_to_i,
                'bm':nr_tm_i,
}

DSargs.varspecs = { 'l':'-k1_on*(1e6*ac/(Na*vims))*l + k1_off*al- k2_on*(1e6*bm/(Na*vmito))*l + k2_off*lb - k5_on*(1e6*bc/(Na*vims))*l + k5_off*bl - k6_on*(1e6*am/(Na*vmito))*l + k6_off*la',
                    'al':'k1_on*(1e6*ac/(Na*vims))*l - k1_off*al - k2_on*(1e6*bm/(Na*vmito))*al + k2_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*ala - k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*alap',
                    'la':'k6_on*(1e6*am/(Na*vmito))*l - k6_off*la  - k5_on*(1e6*bc/(Na*vims))*la + k5_off*bla - k1_on*fa*(1e6*ac/(Na*vims))*la + k1_off*ala - k1_on*fa*(1e6*ac/(Na*vims))*la + k1_off*alap',
                    'lb':'k2_on*(1e6*bm/(Na*vmito))*l - k2_off*lb - k1_on*(1e6*ac/(Na*vims))*lb + k1_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k5_on*fa*(1e6*bc/(Na*vims))*lb + k5_off*blb - k5_on*fa*(1e6*bc/(Na*vims))*lb + k5_off*blbp',
                    'bl':'k5_on*(1e6*bc/(Na*vims))*l - k5_off*bl- k6_on*(1e6*am/(Na*vmito))*bl + k6_off*bla - k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blb- k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blbp',
                    'bla':'k6_on*(1e6*am/(Na*vmito))*bl - k6_off*bla - k8*bla + k7*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) + k5_on*(1e6*bc/(Na*vims))*la - k5_off*bla',
                    'ala':'k1_on*fa*(1e6*ac/(Na*vims))*la - k1_off*ala + k6_on*fa*(1e6*am/(Na*vmito))*al - k6_off*ala -k10*ala + k10*alap',
                    'alap':'k1_on*fa*(1e6*ac/(Na*vims))*la - k1_off*alap + k6_on*fa*(1e6*am/(Na*vmito))*al - k6_off*alap - k10*alap + k10*ala',
                    'blb':'k2_on*fa*(1e6*bm/(Na*vmito))*bl - k2_off*blb + k5_on*fa*(1e6*bc/(Na*vims))*lb - k5_off*blb + k9*blbp - k9*blb',
                    'blbp':'k2_on*fa*(1e6*bm/(Na*vmito))*bl - k2_off*blbp + k5_on*fa*(1e6*bc/(Na*vims))*lb - k5_off*blbp -k9*blbp +k9*blb'}

# initial conditions
DSargs.ics  = {'l':no1_ant,
               'al':0,
               'lb':0,
               'bl':0,
               'bla':0,
               'la':0,
               'ala':0,
               'blb':0,
               'alap':0,
                'blbp':0}

DSargs.tdomain = [0,1]

#ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)    # an instance of the 'Generator' class.
ode  = PyDSTool.Generator.Radau_ODEsystem(DSargs)
traj = ode.compute('polarization')
pd   = traj.sample()



fig = plt.figure(1)
fig.subplots_adjust(right=0.95, left = 0.28, bottom =0.2, top = 0.96)
color = ['0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0']

for i,dmi in enumerate(np.arange(1,10, 0.05)):
        cdm_i = dmi
        ctm_i = 15 - dmi #mM
        print(i,dmi)

        nr_dm_i = 0.45*0.8*6.02*cdm_i*vmito*1e20 #number of molecules Dm
        nr_tm_i = 0.05*6.02*ctm_i*vmito*1e20

        ode.set( pars = {'am':nr_dm_i,
                        'bm':nr_tm_i})


# # #                  # Initial condition
        tmp = ode.compute('pol%3i' % i).sample()    # or specify dt option to sample to sub-sample
        plt.plot(ctm_i/cdm_i, 1e3*k7*(no1_ant - tmp['l'][-1]- tmp['al'][-1]- tmp['bl'][-1]- tmp['la'][-1]- tmp['lb'][-1]- tmp['ala'][-1]- tmp['alap'][-1]- tmp['blbp'][-1]- tmp['blb'][-1]-tmp['bla'][-1])/(Na*vmito) - 1e3*k8*tmp['bla'][-1]/(Na*vmito) ,'.', c=color[0],markersize=1)

plt.ylabel(r'ADP$\rm _i$ - ATP$\rm _m$ flux (mM/s)',labelpad =0.5)
plt.xlabel(r'ATP$\rm _m$/ADP$\rm _m$ ratio ',labelpad =0.5)
plt.ylim(0.15,0.7)
#plt.savefig('my_ant_flux_net_vs_atp_adp_ratiom_last.png',dpi =600)

plt.show()
