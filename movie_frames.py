#from IPython.core.display import display, HTML

#display(HTML("<style>.container {width:100% !important;}</style>"))

#Usage of this script via e.g.
# python3 movie_frames.py "8me-atoms--H1e-14-FUV1e+3-kzz1e-10-1cool-1e4-autogen_singlefluidVIDEO" 0 300 "_singlefluid"
# 
#@args:
#      1: movie_frames.py
#      2: string, simulation files in hard coded folder (runs_2023_co) to make into frames
#      3,4: min input, max output
#      5: additional framelabel
#
#Create movie with
# ffmpeg -framerate 10 -pattern_type glob -i 'frame_*.png' -c:v libx264 -r 10 -pix_fmt yuv420p moviename.mp4
# with autocapping of edge pixels for formatting
# ffmpeg -framerate 10 -pattern_type glob -i 'movieframes/frame_*.png' -c:v libx264 -r 10 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" movieframes/waterworlds2_H1e-2.mp4

import numpy as np
import copy
from scipy import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
#from scipy.special import lambertw
#import scipy.special as sc
from scipy import integrate
from scipy.special import lambertw
matplotlib.get_backend()
#from pyhdf.SD import SD
#from pyhdf.SD import SDC
#from h5py.SD import SD
from decimal import Decimal
import os
import sys

#import pygraphviz as pgv
#from pyflowchart import *
#import pydot
import copy

#from astropy.modeling.models import BlackBody as BB
#from astropy import units as u
#from astropy.visualization import quantity_support
from scipy.integrate import quad

import warnings
warnings.filterwarnings('ignore')

font = {'family': 'serif',
        'color':  'red',
        'weight': 'normal',
        'size': 19,
        }


left1, bottom1, width1, height1 = [0.095, 0.545, 0.88, 0.43]
left2, bottom2, width2, height2 = [0.095, 0.05, 0.88, 0.43]

kb = 1.38e-16
amu= 1.66e-24
G  = 6.678e-8
K_to_eV    = 8.621738e-5;
ev_to_K    = 1./K_to_eV;
sigma_rad = 5.67e-5

rearth   = 6370e5
year      = 365*24*3600

rjupiter = 70000*1e5
rsolar   = 10.*rjupiter
au       = 149.*1e6*1e5

lsolar = 3.83e33 #erg/s
matmo  = 6e21
msun   = 2e33

pi = 3.141592
rpl = 1.25e9
mearth   = 5.98e27
rearth   = 6370e5
micron   = 1e-4
angstr   = 1e-4*micron
debug = 0

lines_O = [
    [6300*angstr, 1.15685197E-14, 22830.7, 8.61213387e+05]]

lines_Op = [
    [834*angstr, 5.78622E-4, 172421.6, 1.32174E15],
[2741*angstr, 3.81198E-13, 58225.3, 4.48777E7],
[3727*angstr, 4.29901E-16, 38575.0, 5.36461E3],
[7320*angstr, 3.76929E-13, 53063.6, 3.11018E7]]

lines_Opp = [
[52*micron, 3.13852E-18, 277.682, 2.5493E3],
[5000*angstr, 3.386678E-14 ,28728.6, 9.66741E5],
[166*angstr, 6.59979E-10, 86632.4, 1.475889E10],
[83.5*angstr, 1.75205E3, 172569.7, 5.405937E21]]

lines_C = [
    [0, 0, 0, 1]]

lines_Cp = [
    [157*micron, 1.78292054E-20, 91.2, 1.38779668E+01],
    [2326*angstr, 1.21471023E-10, 61853.9, 1.20977633e+09],
    [1334*angstr, 2.41304508E-03, 107718.1, 3.74008770E+15]]

lines_Cpp = [
    [1910*angstr, 3.84223024E-10, 75460.8, 1.31478953E+9],
    [977*angstr, 1.79050834E-03, 147263.9, 7.17164800E+14]]

def analytic_wind_solution(mplanet, mparticle, T, radii, gamma_ad = 1.):
    
    mpl= mplanet*6e27
    cs     = np.sqrt(gamma_ad * kb*T/(mparticle*amu))
    rsonic = 0.5*G*mpl/cs**2.
    rrc           = radii/rsonic;
    D             = rrc**(-4.) * exp(4.*(1.-1./rrc)-1.);
    lamb0  = cs/(gamma_ad)*np.sqrt( -lambertw(-D, 0)).real
    lambm1 = cs/(gamma_ad)*np.sqrt( -lambertw(-D, -1)).real
    
    #print("in analytic, rs/rp = ")
    #print(rsonic/rpl)
    
    u_analytic = [ lamb0[i] if radii[i]<rsonic else lambm1[i] for i, rrc in enumerate(rrc)]
    
    return rsonic, cs, u_analytic

def get_numerical_wind_params(data):
    veldcs = np.where(data[:,11]/data[:,15] > 1.)
    
    temper_rs   = data[veldcs,12][0][0]
    rs          = data[veldcs, 0][0][0]
    vel_rs      = data[veldcs, 11][0][0]
    dens_rs     = data[veldcs, 1][0][0]

    return rs, vel_rs, temper_rs, dens_rs

def HOnly_cooling(nX0, nXp, Te) :

    T_HI = 157807.;
    #if (Te < 220):
    #        Te = 220.;
    
    x = 2 * T_HI / Te;
    cooling = 0;

    # Recombination Cooling:
    term = 0 ;
    #if (x < 1e5) :
    #    term = 3.435e-30 * Te * x*x / pow(1 + pow(x / 2.250, 0.376), 3.720);
    #else:
    #    term = 7.562e-27 * pow(Te, 0.42872) ;
        
    term = 3.435e-30 * Te * x*x / pow(1 + pow(x / 2.250, 0.376), 3.720);
    
    cooling += 1.* nXp * term;
    
    # Collisional ionization cooling:
    #term = kb * T_HI * H_collisional_ionization(Te);
    #cooling += 1. *nX[0] * term;

    # HI Line cooling (Lyman alpha):
    term = 7.5e-19 * exp(-0.75 * T_HI / Te); #MC2009
    cooling += 1. * nX0 * term;

    # Free-Free:
    term = 1.426e-27 * 1.3 * sqrt(Te) ;
    cooling += 1. * nXp * term  ;

    return 1.*cooling;


def C_cooling(Te):
    term = 1e-24+3.1e-20*np.exp(-15162/Te)*(1.+pow(Te/2e4, 1.5));
    return term;

def C_cooling2(Te, ne):
    term = 0.
    for ln in lines_C:
        term += ln[1]*np.exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne))
    return term;

def CP_cooling(Te):
    term = 1.5e-23+3.1e-20*np.exp(-45162/Te)*(1.+pow(Te/0.75e4,1.5));
    return term;

def CP_cooling2(Te,ne):
    term = 0.
    for ln in lines_Cp:
        term += ln[1]*np.exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne))
    return term;

def CPP_cooling2(Te,ne):
    term = 0.
    for ln in lines_Cpp:
        term += ln[1]*np.exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne))
    return term;

def O_cooling(Te):
    
    term = 5.5e-24+1.1e-20*np.exp(-30162/Te)*(1.+pow(Te/0.75e4, 0.5));
    
    return term;

def O_cooling2(Te,ne):
    term = 0.
    for ln in lines_O:
        term += ln[1]*np.exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne))
    return term;

def OP_cooling(Te):
    
    term = 5.1e-20*np.exp(-35162/Te)*(1.+pow(Te/0.75e4, 0.5));
    
    return term;

def OP_cooling2(Te,ne):
    term = 0.
    for ln in lines_Op:
        term += ln[1]*np.exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne))
    return term;

def OPP_cooling2(Te,ne):
    term = 0.
    for ln in lines_Opp:
        term += ln[1]*np.exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne))
    return term;



atomicsim=1
imin=0
imax=10
simsim = ""
framelabel = ""
print("argc = " + repr(len(sys.argv)))
for i, arg in enumerate(sys.argv):
    print(f"Argument {i:>6}: {arg}")

if len(sys.argv)>=3:
	simsim  = str(sys.argv[1])
	imin = int(sys.argv[2])
	imax = int(sys.argv[3])

if len(sys.argv)>=5:
	framelabel = str(sys.argv[4])
	print("found framelabel = " + repr(framelabel))
#
#
# Movieframes 1
#
#

rjup = 69799e5
ljup = 1e26
rp= 1189228191.58938 #1.25e9 #0.1*rjup

masses = []
names  = []
rlnames = []
cols = []
spcdict = {'H0': 0, 'S1': 1} #Dummy, properly defined in the code section below	
#masses = [1., 1., 5e-4, 4., 12., 12., 12., 16., 16., 16., 28.]
#masses = [1., 1., 5e-2, 4., 12., 12., 12., 16., 16., 16., 4.]
#masses = [1., 1., 5e-4, 4., 4., 16., 16., 16., 18., 2., 28., 12., 12., 12.]
#names  = ["H0","p+", "e-", "He", "C0", "C+", "C++", "O0", "O+", "O++"]
#names  = ["H0","p+", "eh2"]
#names  = ["S0","S1","S2","S3", "S4","S5","S6", "S7", "S8","S9", "S10"]
#names    = ["S0","S1", "e-", "O0", "Op", "H2O", "H2Op", "H2", "H2p", "O2", "O2p", "OH"]
#rlnames  = ["H ","p+", "e-", "O", "O+",  "H2O", "H2O+", "H2", "H2+", "O2", "O2+", "OH"]
#names  = ["S0","S1", "e-", "He", "C0", "Cp", "Cpp", "O0", "Op", "Opp", "He+"]
#rlnames= ["H0","p+", "e-", "He", "C0", "C+", "C++", "O0", "O+", "O++", "He+"]
#rlnames= ["H0","p+", "e-", "He", "C0", "C+", "C++", "O0", "O+", "O++", "He+"]
#cols   = ["black","red","blue","green","orange","magenta","olive","purple","grey","cyan","magenta"]
#cols   = ["black","black","grey","orange","orange","red","red","red","blue","blue","cyan","cyan", "green", "cyan"]

if atomicsim == 1:
	masses = [1., 1., 5e-4,     4.,4., 12., 12., 12., 16., 16., 16.]
	names  = ["S0","S1", "e-", "He","He+", "C0", "Cp", "Cpp", "O0", "Op", "Opp"]
	rlnames= ["H0","p+", "e-", "He","He+", "C0", "C+", "C++", "O0", "O+", "O++"]
	cols   = ["black","black","dodgerblue","orange","orange","green","green","green","magenta","magenta","magenta","orange","magenta"]
	styles = ['-','--',':',                '-','--',  '-','--',':',           '-','--',':',                 '--','-','-','-', '-', '--', ':', '--']
	spcdict = {'H0': 0, 'p+': 1, 'e-': 2,'He': 3, 'He+': 4,'C0': 5,'Cp': 6,'Cpp': 7,'O0': 8,'Op': 9,'Opp': 10}
else:
	#Water H2 simulations
	masses = [1., 1., 5e-4,     4.,4.,       16., 16., 16.,    18.,2.,28.,           16., 16., 16.]
	names    = ["S0","S1", "e-","He","He+", "O0", "Op", "Opp", "H2O", "H2", "CO", "C0", "Cp", "Cpp"]
	rlnames  = ["H ",r"$\rm p^{+}$", r"$\rm e^{-}$","He",r"$\rm He^{+}$", "O", r"$\rm O^{+}$", r"$\rm O^{++}$",  r"$\rm H_2O$", r"$\rm H_2$", "CO", "C", r"$\rm C^{+}$", r"$\rm C^{++}$"]
	cols   = ["blue","blue","grey","black","black","magenta","magenta","magenta","orange","red","cyan","green", "green", "green"]
	styles = ['-','--',':','-','--','-','--',':','-','-','-','-', '-', '--', ':', '--']
	spcdict = {'H0': 0, 'S1': 1, 'Cpp': 11}
    
times = np.arange(imin,imax,1)
print("times = " + repr(times))
    
#dirr = "../from_william/"
#name = "10me_corepowered5_3s2b_r32_heavyelectrons3"

#dirr = "../c2ray/"
#name = "gj436b_4e6_debug5_3s1b"

#dirr = "runs_co1e-3_4b_fixedbound/"
#name = "8me_H1e-1_UV1e4_lowCO_u1e+2_heavyelectrons10"
#name = "8me_H1e-1_UV1e4_lowCO_u1e-6"
#dirr = "runs_co1e-3_4b_fixedbound/8me_H1e-1_UV1e4_lowCO_u1e-6"
dirr = "runs_2023_co/"
#name = "4me_water_H1e-10_kzz_cocool_movieframes"
#        4me_water_H1e-10_kzz_cocool_movieframes
name= simsim #"8me-atoms--H1e-14-FUV1e+3-kzz1e-10-1cool-1e4-autogen_singlefluidVIDEO"
#name="8me-atoms--H1e-14-FUV1e+3-kzz1e-10-1cool-1e4-autogen_electronkzzpotentialallowed"

file_length = 699 #260 #170
    
for i,time in enumerate(times):    
    #print("starting frame "+repr(time))
    fig, axs = plt.subplots(2, 3, figsize=(24, 16))
    ax01a = axs[0,1].twinx()
    inset = axs[0,1].inset_axes([0.45,0.1,0.25,0.40])
    
    for spl in axs.flat:
        spl.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
        spl.tick_params(axis='x', which='major', labelsize=25, size=10,direction="in");
        spl.tick_params(axis='y', which='major', labelsize=28, size=10,direction="in");
        spl.yaxis.set_ticks_position('both')
        spl.xaxis.set_ticks_position('both') 
        spl.set_xlabel(r'$\rm Radius \; [\rm R_P]$', fontsize=17, labelpad=-15)
        spl.set_xlim([0.95, 2.05])

    axs[0,0].set_xlabel(r'$\rm Radius \; [\rm R_P]$', fontsize=16, labelpad=0)
    axs[0,0].tick_params(axis='x', which='major', labelsize=20);
    axs[0,1].yaxis.set_ticks_position('left')
    axs[0,2].yaxis.set_ticks_position('right')
    inset.yaxis.set_ticks_position('right')
    for spl in [ax01a, inset]:
        spl.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
    #ax01a.tick_params(axis='x', which='major', labelsize=25, size=10,direction="in");
        spl.tick_params(axis='y', which='major', labelsize=23, size=10,direction="in");
    #ax01a.yaxis.set_ticks_position('both')
    #ax01a.xaxis.set_ticks_position('both')
    #axs[0,2].ax01a.tick_params(axis='y', which='major', labelsize=23, size=10,direction="in");
    spl.tick_params(axis='both', which='major', labelsize=20, size=10,direction="in");

    axs[0,0].set_ylim([1e-5, 1e20])    
    axs[0,0].set_xlim([0.90, 2.05])    
    axs[0,2].set_ylim([1e1, 1e11])
    axs[0,1].set_ylim([1e2, 3e4])
    ax01a.set_ylim([1e1, 1.5e7])
    axs[1,0].set_ylim([1e-30, 1e0])
    axs[1,1].set_ylim([1e-30, 1e0])
    axs[1,2].set_ylim([0.9e-14, 1.5e3])
    inset.set_xlim([9.,10.])
    inset.set_ylim([1e5,3.1e6])

    axs[0,0].set_title(r"$\rm Number\,density\;[cm^{-3}]$",fontsize=18)
    axs[0,2].set_title(r"$\rm Mass\,flux\;[g\;s^{-1}]$",fontsize=18)
    axs[0,1].set_title(r"$\rm Temperature\,[K](grey,left)\; /\; Velocity\,[cm\;s^{-1}](right)$",fontsize=17)
    
    axs[1,0].set_title(r"$\rm Heating\,[erg\;s^{-1}\;cm^{-3}]$",fontsize=18)
    axs[1,1].set_title(r"$\rm Cooling\,[erg\;s^{-1}\;cm^{-3}]$",fontsize=18)
    axs[1,2].set_title(r"$\rm Total\,pressure\;[bar]$",fontsize=18)
    #ax01a.set_title(r"$\rm Velocity\,[cm/s]$",fontsize=18)
    
    numdensies = []

    dat = np.loadtxt(dirr+"output_"+name+"_"+names[0] +"_t"+ repr(time)+".dat",  skiprows=2, max_rows=file_length-3)
    totalpressure = dat[:,10] * 0.

    for j,spc in enumerate(names):
        
        dat = np.loadtxt(dirr+"output_"+name+"_"+spc +"_t"+ repr(time)+".dat",  skiprows=2, max_rows=file_length-3)
        #diag = np.loadtxt(dirr+"diagnostic_"+name+"_t"+ repr(time)+".dat",  skiprows=2, max_rows=file_length-3)
        
        r    = dat[:,0]/rp  #years
        dens = dat[:,1]
        totalpressure += dat[:,10]
        numdensies.append(dens/(masses[j]*amu))
        #rsonic_ana, cs_ana, u_ana = analytic_wind_solution(8., 1., 10000., dat[:,0], 1.0)
	
        axs[0, 0].semilogy(r, dens/(masses[j]*amu), lw=3., label=rlnames[j], c=cols[j], ls=styles[j])
        axs[0, 2].loglog(r, 12.566*dat[:,2]*r*r*rp*rp, lw=3.,      label=rlnames[j], c=cols[j], ls=styles[j])
        #axs[0, 1].loglog(r, dat[:,12], lw=3.,      label=rlnames[j], c=cols[j], ls=styles[j]) #New temperature curve is grey
        axs[0, 1].loglog(r, dat[:,12], lw=2.,      label=rlnames[j], c='grey', ls=styles[j])
        
        axs[1, 0].loglog(r, dat[:,22], lw=3.,      label=rlnames[j], c=cols[j], ls=styles[j])
        #axs[1, 1].loglog(r, dat[:,21], lw=3.,      label=rlnames[j], c=cols[j], ls=styles[j])
        if 0==0:
           ax01a.loglog(r, dat[:,11], lw=3.,      label=rlnames[j], c=cols[j], ls=styles[j])
           inset.semilogy(r, dat[:,11], lw=3.,    label=rlnames[j], c=cols[j], ls=styles[j])
           #inset.semilogy(r, np.log10(dat[:,11]), lw=3.,    label=rlnames[j], c='grey', ls=styles[j])
           ax01a.loglog(r, -dat[:,11], lw=1.,      label="", c=cols[j], ls=styles[j])
        else:
           ax01a.loglog(r, dat[:,11], lw=3.,      label=rlnames[j], c=cols[j], ls='--')
           inset.plot(r, dat[:,11], lw=3.,      label=rlnames[j], c=cols[j], ls='--')
           ax01a.loglog(r, -dat[:,11], lw=1.,      label="", c=cols[j], ls='--')
        
        axs[0,0].text(x=1.6,y=1e14,s='Frame='+repr(time),fontsize=18)

    axs[1,2].loglog(r, totalpressure/1e6, lw=3., c='k')

    nH = numdensies[spcdict['H0']]
    nP = numdensies[spcdict['p+']]
    nE = numdensies[spcdict['e-']] 
    nC = numdensies[spcdict['C0']]
    nCp= numdensies[spcdict['Cp']]
    nCpp= numdensies[spcdict['Cpp']]
    nO = numdensies[spcdict['O0']] 
    nOp= numdensies[spcdict['Op']]
    nOpp= numdensies[spcdict['Opp']]
    dat = np.loadtxt(dirr+"output_"+name+"_e-_t"+ repr(time)+".dat",  skiprows=2, max_rows=file_length-3)
    Te = dat[:,12] 

    coolH  = nE  * HOnly_cooling(nH,nP, Te)
    coolC  = nC  * nE * C_cooling2(Te,nE); 
    coolCp = nCp * nE * CP_cooling2(Te,nE); 
    coolCpp = nCpp * nE * CPP_cooling2(Te,nE); 
    coolO  = nO  * nE * O_cooling2(Te,nE); 
    coolOp = nOp * nE * OP_cooling2(Te,nE); 
    coolOpp = nOpp * nE * OPP_cooling2(Te,nE); 
    nullc = 0. * CPP_cooling2(Te,nE); 

    cools=[]
    if atomicsim == 1:
        cools = [coolH, nullc,nullc,nullc,nullc,coolC, coolCp, coolCpp, coolO, coolOp, coolOpp,  nullc, nullc, nullc, nullc, nullc, nullc, nullc]
    else:
        cools = [coolH, nullc,nullc,nullc,nullc, coolO, coolOp, coolOpp, nullc, nullc, nullc,coolC,coolCp,coolCpp, nullc, nullc, nullc, nullc]
    
    for j,spc in enumerate(names):
        axs[1, 1].loglog(r, cools[j], lw=3.,      label=rlnames[j], c=cols[j], ls=styles[j])

    for ax in [axs[0,0]]:
        ax.legend(fontsize=14, loc='upper left', ncol=4)
    
    for ax in [axs[1,0], axs[1,1] ]:
        ax.legend(fontsize=14)

    for ax in [axs[0,2]]:
        ax.legend(fontsize=14, loc='upper right', title='', bbox_to_anchor=(1., 1.))

    if time>999:
        zeros = ""
    else:
        if time>99:
           zeros = "0"
        else:
            if time>9:
               zeros = "00"
            else:
               zeros = "000"
            
    print("saving frame "+repr(time))
    
    fig.savefig("movieframes/frame"+framelabel+"_"+zeros+repr(time)+".png",bbox_inches='tight')
