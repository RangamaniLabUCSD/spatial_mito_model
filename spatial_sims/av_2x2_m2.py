import PyDSTool
import matplotlib.pyplot as plt
import pickle
import numpy as np
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
from scipy import stats
from matplotlib.lines import Line2D
import matplotlib as mpl

'''
This script runs the ode model, and plot the results along the average traces of the spatial model.
Fig 8B in the paper
08/11/23
GCG
'''
params = {'axes.labelsize': 9,
           'axes.titlesize': 9,
          'legend.fontsize': 9,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
            'figure.figsize': (2.5,2.5)}
mpl.rcParams.update(params)
#

Na = 6.02214e23

#Mito Volume
vmito = 0.010e-15 #IM volume  ###########
vims = 0.019e-15 #IMS volume ########
vcube = 0.061e-15 #volume of the cytosol #########

cdm_i = 2 #10mM
ctm_i = 13 #mM
cdo_i = 0.1 #10mM
cto_i = 6.5 #mM

nr_dm_m = np.round(0.45*0.8*6.02*cdm_i*vmito*1e20,decimals =0) #number of molecules Dm
nr_tm_m = np.round(0.05*6.02*ctm_i*vmito*1e20,decimals =0)
#print(dm_i,tm_i)
nr_do_i = np.round(0.45*6.02*cdo_i*vims*1e20,decimals =0)
nr_to_i = np.round(0.05*6.02*cto_i*vims*1e20,decimals =0)
nr_to_c = np.round(0.05*6.02*cto_i*vcube*1e20,decimals =0)

print('matrix',nr_dm_m,nr_tm_m)
print('ims',nr_do_i,nr_to_i,nr_to_c,0.05*6.5, 0.1*0.45)

#Number of proteins
no1_ant = 10718 #number of ANTS #############
no_atp = 150 #number of atpases #######
n_porin = 5985 #### Number of VDACS #########

k7 = 92 #0.74
k8 = 3.5 #0.05s
k9 = 5.8 #0.37
k10 = 4.8#0.47
cm = 3e3
a12 = 100
a21 = 40.0
a23 = 5.0
a32 = 5e3
ac = nr_do_i #clamp*Na*vims
kp = 1.0
DSargs = PyDSTool.args(name='ex')

DSargs.pars = { 'k1_on':'4.0',#          #KDo
                'k1_off':'100', #         #KDo
                'k2_on': '6.4',#         #KTi
		'k2_off':'4e4',#     #KTi
		'k5_on':'0.4',           #KTo
		'k5_off':'200.0',           #KTo
		'k6_on':'4.0',           #KDi
		'k6_off':'4e4',      #KDi
		'k7':k7,           #kp
		'k8':k8,              #kcp
		'k9':k9,             #kt
		'k10':k10,            #kd
		'no_ant':no1_ant ,
		'a12':a12, #s-1
		'a21':a21,#s-1
		'a65':'3969',           #s-1
		'a56':'1e3',             #s-1,
		'a16':'146516',           #s-1
		'a61':'33989',           #s-1
		'a25':'0',     #s-1
		'a52':'0.0',          #s-1
		'a54':'1e2',          #s-1
		'a45':'1e2',          #s-1
		'a43':'2.0',         #uM-1s-1
		'a34':'1e2',
		'a32':a32,          #s-1
		'a23':a23,         #uM-1s-1
		'no1_atp': no_atp, #number of atphase
		'vims':vims,
		'vcube':vcube,
		'vmito':vmito,
		'kp':kp,
		'n_porin':n_porin,
		'Na':'6.02214e23',
		'fa':'0.5',
		'ac':nr_do_i}
#DSargs.fnspecs  = {'alb': (['l','al','lb','bl','bla','la','ala','blb','alap','blbp'], 'no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp') }

DSargs.varspecs = { 'am':'-k6_on*(1e6*am/(Na*vmito))*bl + k6_off*bla - k6_on*(1e6*am/(Na*vmito))*l + k6_off*la - 2.0*k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*ala + k6_off*alap + a34*h3es - a43*(1e6*am/(Na*vmito))*h3eo',
                    'bm':'-k2_on*(1e6*bm/(Na*vmito))*al + k2_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k2_on*(1e6*bm/(Na*vmito))*l + k2_off*lb - 2.0*k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blb + k2_off*blbp -a23*h3e*(1e6*bm/(Na*vmito)) + a32*h3es',
                    'bc':'-k5_on*(1e6*bc/(Na*vims))*l + k5_off*bl - k5_on*(1e6*bc/(Na*vims))*la + k5_off*bla- 2.0*k5_on*fa*(1e6*bc/(Na*vims))*lb + k5_off*blb + k5_off*blbp - kp*(1e6*bc/(Na*vims))*n_porin + kp*(1e6*bo/(Na*vcube))*n_porin',
                    'bo':' n_porin*kp*(1e6*bc/(Na*vims)) - n_porin*kp*(1e6*bo/(Na*vcube))',
                    'l':'-k1_on*(1e6*ac/(Na*vims))*l + k1_off*al- k2_on*(1e6*bm/(Na*vmito))*l + k2_off*lb - k5_on*(1e6*bc/(Na*vims))*l + k5_off*bl - k6_on*(1e6*am/(Na*vmito))*l + k6_off*la',
                    'al':'k1_on*(1e6*ac/(Na*vims))*l - k1_off*al - k2_on*(1e6*bm/(Na*vmito))*al + k2_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*ala - k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*alap',
                    'la':'k6_on*(1e6*am/(Na*vmito))*l - k6_off*la  - k5_on*(1e6*bc/(Na*vims))*la + k5_off*bla - k1_on*fa*(1e6*ac/(Na*vims))*la + k1_off*ala - k1_on*fa*(1e6*ac/(Na*vims))*la + k1_off*alap',
                    'lb':'k2_on*(1e6*bm/(Na*vmito))*l - k2_off*lb - k1_on*(1e6*ac/(Na*vims))*lb + k1_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k5_on*fa*(1e6*bc/(Na*vims))*lb + k5_off*blb - k5_on*fa*(1e6*bc/(Na*vims))*lb + k5_off*blbp',
                    'bl':'k5_on*(1e6*bc/(Na*vims))*l - k5_off*bl- k6_on*(1e6*am/(Na*vmito))*bl + k6_off*bla - k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blb- k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blbp',
                    'bla':'k6_on*(1e6*am/(Na*vmito))*bl - k6_off*bla - k8*bla + k7*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) + k5_on*(1e6*bc/(Na*vims))*la - k5_off*bla',
                     'ala':'k1_on*fa*(1e6*ac/(Na*vims))*la - k1_off*ala + k6_on*fa*(1e6*am/(Na*vmito))*al - k6_off*ala -k10*ala + k10*alap',
                     'alap':'k1_on*fa*(1e6*ac/(Na*vims))*la - k1_off*alap + k6_on*fa*(1e6*am/(Na*vmito))*al - k6_off*alap - k10*alap + k10*ala',
                     'blb':'k2_on*fa*(1e6*bm/(Na*vmito))*bl - k2_off*blb + k5_on*fa*(1e6*bc/(Na*vims))*lb - k5_off*blb + k9*blbp - k9*blb',
                     'blbp':'k2_on*fa*(1e6*bm/(Na*vmito))*bl - k2_off*blbp + k5_on*fa*(1e6*bc/(Na*vims))*lb - k5_off*blbp -k9*blbp +k9*blb',
                     'eo':'-a65*eo + a56*(no1_atp-eo-ei-h3eo-h3es-h3e)+ a16*ei - a61*eo',
                     'ei':'-a16*ei + a61*eo - a12*ei + a21*h3e',
                     'h3eo':'-a45*h3eo + a54*(no1_atp-eo-ei-h3eo-h3es-h3e) + a34*h3es - a43*h3eo*(1e6*am/(Na*vmito))',
                   'h3es':'a43*(1e6*am/(Na*vmito))*h3eo - a34*h3es + a23*h3e*(1e6*(bm)/(Na*vmito)) - a32*h3es',
                   'h3e':'-a23*h3e*(1e6*(bm)/(Na*vmito)) + a32*h3es - a25*h3e + a52*(no1_atp-eo-ei-h3eo-h3es-h3e) + a12*ei - a21*h3e'}

# initial conditions
DSargs.ics  = {'am':nr_dm_m,
	       'bm':nr_tm_m,
	       'bc':nr_to_i,
	       'bo':nr_to_c,
	       'l':no1_ant,
                'al':0,
                'lb':0,
                'bl':0,
                'bla':0,
                'la':0,
                'ala':0,
                'blb':0,
                'alap':0,
                'blbp':0,
                'eo':no_atp,
                'ei':0,
                'h3eo':0,
                'h3es':0,
                'h3e':0}

DSargs.tdomain = [0,0.3]
#DSargs.algparams   =   {'init_step':0.00001, 'strictdt':True , 'max_pts': 5000000}

#ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)    # an instance of the 'Generator' class.
ode  = PyDSTool.Generator.Radau_ODEsystem(DSargs)
traj = ode.compute('polarization')
#print '  ... finished in %.3f seconds.\n' % (clock()-start)
pd   = traj.sample()


time = np.genfromtxt("./mito2/high_curva/time_10c_1e8_m2", dtype = float)
av_var = pickle.load(open('./mito2/high_curva/av_10c_1e8_m2','rb'))

fig = plt.figure(1)
fig.subplots_adjust(right=0.97, left = 0.16, bottom =0.15, top = 0.90, wspace=0.7, hspace = 0.2)
ax = plt.subplot(2,2,1)

ax1 = plt.subplot(2,2,1)
plt.plot(1e3*time[:],av_var[1,:],'b')
plt.plot(1e3*pd['t'], pd['am'],c='k')
ax1.locator_params(axis='y',nbins=2)
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel(r'#ADP$_{\rm matrix}$')
ax1.locator_params(axis='x',nbins=8)
ax1.xaxis.set_major_formatter(plt.NullFormatter())
ax1.set_xticks([0,100,200,300])

ax2 = plt.subplot(2,2,2)
plt.plot(1e3*time[:],av_var[2,:],'b')
plt.plot(1e3*pd['t'], pd['bm'],c='k')
ax2.locator_params(axis='y',nbins=3)
ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel(r'#ATP$_{\rm matrix}$')
ax2.locator_params(axis='x',nbins=8)
ax2.xaxis.set_major_formatter(plt.NullFormatter())
ax2.set_xticks([0,100,200,300])

ax3 = plt.subplot(2,2,3)
plt.plot(1e3*time[:],av_var[3,:],'b')
plt.plot(1e3*pd['t'], pd['bc'],c='k')
ax3.locator_params(axis='y',nbins=3)
ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel(r'#ATP$_{\rm IMS}$')
ax3.locator_params(axis='x',nbins=8)
plt.xlabel('Time (msec)')
ax3.set_xticks([0,100,200,300])

ax4 = plt.subplot(2,2,4)
plt.plot(1e3*time[:],av_var[29,:],'b')
plt.plot(1e3*pd['t'],pd['bo'],'k')
ax4.locator_params(axis='y',nbins=3)
ax4.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Time (msec)')
plt.ylabel(r'# ATP$_{\rm cytosol}$')
ax4.locator_params(axis='x',nbins=8)
ax4.locator_params(axis='x',nbins=8)
ax4.set_xticks([0,100,200,300])


plt.savefig('av_mito_2.png',dpi=600,bbox_inches='tight')
plt.show()
