import PyDSTool
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
import matplotlib as mpl
from matplotlib.lines import Line2D

'''
Nr of molecules vs time, for fast ants, fig 5 in the paper
08.11.23
GCG
'''

params = {'axes.labelsize': 8,
           'axes.titlesize': 6,
          'legend.fontsize': 6,
           'xtick.labelsize': 6,
           'ytick.labelsize': 6,
            'figure.figsize': (4,4.3)}
mpl.rcParams.update(params)

Na = 6.02214e23
vmito = 0.016e-15 #vol matrix lt
vims = 0.021e-15 # vol IMS
vcube = 0.306e-15 #vol cyto

#---
cdm_i = 2 #mM adp concentration matrix
ctm_i = 13 #mM atp concentration matrix
cdo_i = 0.1 #mM adp concentration cyto/ims 
cto_i = 6.5 #mM atp concentration cyto/ims

nr_dm_m = 0.45*0.8*6.02*cdm_i*vmito*1e20 #number of adp molecules in the matrix
nr_tm_m = 0.05*6.02*ctm_i*vmito*1e20 # number of atp molecules in the matrix 
nr_do_i = 0.45*6.02*cdo_i*vims*1e20 #number of adp molecules in the ims
nr_to_i = 0.05*6.02*cto_i*vims*1e20 #number of atp molecules in the ims
nr_to_c = 0.05*6.02*cto_i*vcube*1e20 #number of atp molecules in the cyto


no1_ant = 16471 #nr of ants
no_atp = 267 #number of atphase
k7 = 9.2 #0.74
k8 = 0.35 #0.05
k9 = 0.58 #0.37
k10 = 0.48#0.47
a12 = 100 #s-1
a21 = 40.0
a23 = 5.0
a32 = 5e3
a65 = 3969
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
		        'a16':'146516',           #s-1
		        'a61':'33989',            #s-1
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
}

DSargs.varspecs = { 'am':'-k6_on*(1e6*am/(Na*vmito))*bl + k6_off*bla - k6_on*(1e6*am/(Na*vmito))*l + k6_off*la - 2.0*k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*ala + k6_off*alap + a34*h3es - a43*(1e6*am/(Na*vmito))*h3eo',
                    'bm':'-k2_on*(1e6*bm/(Na*vmito))*al + k2_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k2_on*(1e6*bm/(Na*vmito))*l + k2_off*lb - 2.0*k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blb + k2_off*blbp -a23*h3e*(1e6*bm/(Na*vmito)) + a32*h3es',
                    'bo':'-kp*(1e6*bo/(Na*vcube))*n_porin + kp*(1e6*bims/(Na*vims))*n_porin',
                    'bims':'-k5_on*(1e6*bims/(Na*vims))*l + k5_off*bl - k5_on*(1e6*bims/(Na*vims))*la + k5_off*bla- 2.0*k5_on*fa*(1e6*bims/(Na*vims))*lb + k5_off*blb + k5_off*blbp - kp*(1e6*bims/(Na*vims))*n_porin + kp*(1e6*bo/(Na*vcube))*n_porin',
                    'l':'-k1_on*(1e6*ac/(Na*vims))*l + k1_off*al- k2_on*(1e6*bm/(Na*vmito))*l + k2_off*lb - k5_on*(1e6*bims/(Na*vims))*l + k5_off*bl - k6_on*(1e6*am/(Na*vmito))*l + k6_off*la',
                    'al':'k1_on*(1e6*ac/(Na*vims))*l - k1_off*al - k2_on*(1e6*bm/(Na*vmito))*al + k2_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*ala - k6_on*fa*(1e6*am/(Na*vmito))*al + k6_off*alap',
                    'la':'k6_on*(1e6*am/(Na*vmito))*l - k6_off*la  - k5_on*(1e6*bims/(Na*vims))*la + k5_off*bla - k1_on*fa*(1e6*ac/(Na*vims))*la + k1_off*ala - k1_on*fa*(1e6*ac/(Na*vims))*la + k1_off*alap',
                    'lb':'k2_on*(1e6*bm/(Na*vmito))*l - k2_off*lb - k1_on*(1e6*ac/(Na*vims))*lb + k1_off*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) - k5_on*fa*(1e6*bims/(Na*vims))*lb + k5_off*blb - k5_on*fa*(1e6*bims/(Na*vims))*lb + k5_off*blbp',
                    'bl':'k5_on*(1e6*bims/(Na*vims))*l - k5_off*bl- k6_on*(1e6*am/(Na*vmito))*bl + k6_off*bla - k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blb- k2_on*fa*(1e6*bm/(Na*vmito))*bl + k2_off*blbp',
                    'bla':'k6_on*(1e6*am/(Na*vmito))*bl - k6_off*bla - k8*bla + k7*(no_ant-l-al-la-lb-bl-bla-ala-blb-alap-blbp) + k5_on*(1e6*bims/(Na*vims))*la - k5_off*bla',
                    'ala':'k1_on*fa*(1e6*ac/(Na*vims))*la - k1_off*ala + k6_on*fa*(1e6*am/(Na*vmito))*al - k6_off*ala -k10*ala + k10*alap',
                    'alap':'k1_on*fa*(1e6*ac/(Na*vims))*la - k1_off*alap + k6_on*fa*(1e6*am/(Na*vmito))*al - k6_off*alap - k10*alap + k10*ala',
                    'blb':'k2_on*fa*(1e6*bm/(Na*vmito))*bl - k2_off*blb + k5_on*fa*(1e6*bims/(Na*vims))*lb - k5_off*blb + k9*blbp - k9*blb',
                    'blbp':'k2_on*fa*(1e6*bm/(Na*vmito))*bl - k2_off*blbp + k5_on*fa*(1e6*bims/(Na*vims))*lb - k5_off*blbp -k9*blbp +k9*blb',
                    'eo':'-a65*eo + a56*(no1_atp-eo-ei-h3eo-h3es-h3e)+ a16*ei - a61*eo',
                    'ei':'-a16*ei + a61*eo - a12*ei + a21*h3e',
                    'h3eo':'-a45*h3eo + a54*(no1_atp-eo-ei-h3eo-h3es-h3e) + a34*h3es - a43*h3eo*(1e6*am/(Na*vmito))',
                    'h3es':'a43*(1e6*am/(Na*vmito))*h3eo - a34*h3es + a23*h3e*(1e6*(bm)/(Na*vmito)) - a32*h3es',
                    'h3e':'-a23*h3e*(1e6*(bm)/(Na*vmito)) + a32*h3es - a25*h3e + a52*(no1_atp-eo-ei-h3eo-h3es-h3e) + a12*ei - a21*h3e'}

# initial conditions
DSargs.ics  = {'am':nr_dm_m,
	           'bm':nr_tm_m,
	           'bims':nr_to_i, #IMS
	           'bo':nr_to_c, #outside
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

DSargs.tdomain = [0,1]

#ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)    # an instance of the 'Generator' class.
ode  = PyDSTool.Generator.Radau_ODEsystem(DSargs)
traj = ode.compute('polarization')
pd   = traj.sample()

fig = plt.figure(1)
fig.subplots_adjust(right=0.95, left = 0.13, bottom =0.1, top = 0.95,wspace=0.5, hspace = 0.6)
color = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','0.93']
for n,j in enumerate(np.arange(1,11,1)):

    k7 = j*9.2 #0.74
    k8 = j*0.35 #0.05
    k9 = j*0.58 #0.37
    k10 = j*0.48#0.47
    #print(n,j)

    ode.set( pars = {'k7':k7,                 #kp
                     'k8':k8,                 #kcp
                     'k9':k9,                 #kt
                     'k10':k10})
# # #                  # Initial condition
    tmp = ode.compute('pol%3i' % n).sample()    # or specify dt option to sample to sub-sample
    the_filename = 'tmp_fast_r_last_'+str(np.round(j,decimals=0))
    #with open(the_filename, 'wb') as f:#
    #    pickle.dump(tmp, f)
    tmp_1 = pickle.load(open('./tmp_fast_r_last_1','rb'))

    ax = plt.subplot(3,2,1)

    ax3 = plt.subplot(3,2,1)
    plt.plot(1e3*tmp['t'], tmp['am'],c=color[n])
    ax3.locator_params(axis='y',nbins=2)
    ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Time (msec)')
    plt.ylabel(r'#ADP$_{\rm matrix}$')
    ax3.locator_params(axis='x',nbins=8)


    ax2 = plt.subplot(3,2,2)
    plt.plot(1e3*tmp['t'], tmp['bm'],c=color[n])
    ax2.locator_params(axis='y',nbins=3)
    ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Time (msec)')
    plt.ylabel(r'#ATP$_{\rm matrix}$')
    ax2.locator_params(axis='x',nbins=8)


    ax1 = plt.subplot(3,2,3)
    plt.plot(1e3*tmp['t'], tmp['bims'],c=color[n])
    ax1.locator_params(axis='y',nbins=3)
    ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylim(3.5e3,4.7e3)
    plt.ylabel(r'#ATP$_{\rm IMS}$')
    ax1.locator_params(axis='x',nbins=8)
    plt.xlabel('Time (msec)')


    ax2 = plt.subplot(3,2,4)
    plt.plot(1e3*tmp['t'],tmp['bo'],color[n])#-tmp['bo'][0]
    ax2.locator_params(axis='y',nbins=3)
    ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Time (msec)')
    plt.ylabel(r'$\Delta$ #ATP$_{\rm cytosol}$')
    ax2.locator_params(axis='x',nbins=8)

    ax6 = plt.subplot(3,2,5)
    plt.plot(j,100*(tmp['bo'][-1]-tmp_1['bo'][-1])/tmp_1['bo'][-1],marker='o',c=color[n])
    plt.xlim(0,11)
    #plt.ylim(-1,6)
    ax6.locator_params(axis='y',nbins=3)
    ax6.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Mutiplicative factor r')
    plt.ylabel("% ATP increment" "\n after 1 sec",labelpad=0.3)
    ax6.locator_params(axis='x',nbins=8)

circ1 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.1",mec='0.1')
circ2 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.2",mec='0.2')
circ3 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.3",mec='0.3')
circ4 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.4",mec='0.4')
circ5 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.5",mec='0.5')
circ6 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.6",mec='0.6')
circ7 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.7",mec='0.7')
circ8 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.8",mec='0.8')
circ9 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.9",mec='0.9')
circ10 = Line2D([0], [0], linestyle="none", marker="o",  markersize=3, markerfacecolor="0.93",mec='0.93')
plt.legend((circ1, circ2,circ3,circ4,circ5,circ6,circ7,circ8,circ9,circ10), ("r = 1","r = 2","r = 3","r = 4","r = 5","r = 6","r = 7","r = 8","r = 9","r = 10"), numpoints=1, loc="upper left", frameon = True,bbox_to_anchor=(1.7, 1.15))

#plt.savefig('param_sweep_ant_speed_atpc_last1.png', format = 'png',dpi=600,bbox_inches='tight')

plt.show()
