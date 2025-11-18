
import numpy as np
import copy
from scipy import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
from scipy.special import lambertw
import scipy.special as sc
from scipy import integrate
matplotlib.get_backend()
from pyhdf.SD import SD
from pyhdf.SD import SDC
#from h5py.SD import SD
from decimal import Decimal
import os
import matplotlib.transforms as mtransforms
import math
from decimal import Decimal

#import pygraphviz as pgv
#from pyflowchart import *
#import pydot
#
#from p_winds import tools, parker, hydrogen, helium, transit, lines

import warnings
warnings.filterwarnings('ignore')

font = {'family': 'serif',
        'color':  'red',
        'weight': 'normal',
        'size': 19,
        }

import logging
logging.getLogger().setLevel(logging.CRITICAL)

left1, bottom1, width1, height1 = [0.095, 0.545, 0.88, 0.43]
left2, bottom2, width2, height2 = [0.095, 0.05, 0.88, 0.43]

kb = 1.38e-16
amu= 1.66e-24
G  = 6.678e-8
K_to_eV    = 8.621738e-5;
ev_to_K    = 1./K_to_eV;
sigma_rad = 5.670374419e-5

pi = 3.141592
rpl = 9.45e9
mearth   = 5.98e27
rearth   = 6370e5
rsolar   = 695510e5;

def parker_massloss(T, mplanet,mparticle, rplanet, dens0, gamma_ad):
    
    rp     = rplanet# * rearth
    mpl    = mplanet * mearth
    cs     = np.sqrt(gamma_ad*kb*T/(mparticle*amu))
    rsonic = 0.5*G*mpl/cs**2.
    
    return 4.*3.141592*rsonic**2.*cs*dens0*exp(-G*mpl/cs**2.*(1/rp-1/rsonic)-0.5)


def cranmer_wind_solution(mplanet, mparticle, T, radii):
    
    mpl= mplanet*6e27
    cs     = np.sqrt(kb*T/(mparticle*amu))
    rsonic = 0.5*G*mpl/cs**2.
    rrc           = radii/rsonic;
    print(" in wind solution, len(radii)  = " + repr(len(radii)))
    print(" in wind solution, len(rrc)  = " + repr(len(rrc)))
    D             = rrc**(-4.) * exp(4.*(1.-1./rrc)-1.);
    lamb0  = cs*np.sqrt( -lambertw(-D, 0)).real
    lambm1 = cs*np.sqrt( -lambertw(-D, -1)).real
    
    print("in analytic, rs/rp = ")
    print(rsonic/rpl)
    
    u_analytic = [ lamb0[i] if radii[i]<rsonic else lambm1[i] for i, rrc in enumerate(rrc)]
    print(" in wind solution, len(u)  = " + repr(len(u_analytic)))
    
    return rsonic, cs, u_analytic


def analytic_wind_solution(mplanet, mparticle, T, radii):
    
    mpl= mplanet*6e27
    cs     = np.sqrt(kb*T/(mparticle*amu))
    rsonic = 0.5*G*mpl/cs**2.
    rrc           = radii/rsonic;
    print(" in wind solution, len(radii)  = " + repr(len(radii)))
    print(" in wind solution, len(rrc)  = " + repr(len(rrc)))
    D             = rrc**(-4.) * np.exp(4.*(1.-1./rrc)-1.);
    lamb0  = cs*np.sqrt( -lambertw(-D, 0)).real
    lambm1 = cs*np.sqrt( -lambertw(-D, -1)).real
    
    print("in analytic, rs/rp = ")
    print(rsonic/rpl)
    
    u_analytic = [ lamb0[i] if radii[i]<rsonic else lambm1[i] for i, rrc in enumerate(rrc)]
    print(" in wind solution, len(u)  = " + repr(len(u_analytic)))
    
    #return rsonic, cs, u_analytic
    return u_analytic
def get_filelength(dirr, file):
    
    count=0
    with open(dirr + file, 'r') as fp:
        for count, line in enumerate(fp):
            pass
    return count-3

# Takes the hydrogen data, FUV flux, planet mass in earth masses and absorption radius of FUV flux. Default cp is max. cp i.e. that of pure hydrogen.
def get_energylimit_massloss_bad(F, mplanet, rabs, data, cp=1.467258e+07, badness=0):
    
    iso_cs     = data[np.where(data[:,11]/data[:,15] > 1.),15]
    iso_T      = data[np.where(data[:,11]/data[:,15] > 1.),12]
    
    rp        = data[1,0]
    cs        = iso_cs[0][0] #data[-1,15]
    cs_low    = data[10,15]
    Thigh     = iso_T[0][0] #data[-1,12]
    Tlow      = data[10,12]
    rs_approx = 0.5*G*mplanet*mearth/cs**2.
    
    if badness == 0:
        return 3.141592*F* (rabs)**2./ (G*mplanet*mearth * (1./rp-1./rs_approx) + 0.5*cs**2 + gamma_ad*cv*(Thigh-Tlow))
    elif badness == 1:
        return 3.141592*F* (rabs)**2./ (G*mplanet*mearth * (1./rp-1./rs_approx))
    elif badness == 2:
        return 3.141592*F* (rabs)**3. / (G*mplanet*mearth)
    elif badness == 3:
        return 3.141592*F* (rabs)**2. * rp / (G*mplanet*mearth)
    else:
        return 3.141592*F* rp**3. / (G*mplanet*mearth)
    
def get_energylimit_massloss_onlyplanetradius(F, mplanet, rplanet, cp=1.467258e+07, badness=0):
    return 3.141592*F* (rplanet)**3. / (G*mplanet*mearth)
    


#Returns ratio of mdots, input are escaping mass flux, particle mass ratio, gravitational parameter and muexp factor to get over the pole
def get_analytic_fractionation(phi, mu, lambda0, muexp=0.): 
    
    nom     = mu*phi - (mu-1.)
    denom   = mu*phi - (mu-1.)*(mu)**muexp*np.exp(-lambda0* (mu*phi - (mu-1.) ) ) #
    
    return nom/denom


def get_analytic_fractionation_hunten(phi, mu, lambda0): #The analytic crossover mass by Hunten
    
    return 1.-(mu-1.)/(mu*phi)


def get_data_fractionation(data_h, data_he ):
    
        m0 = 1.*amu
        m1 = 4.*amu
    
        #data = np.loadtxt(run+"H2_t-1.dat", skiprows=2, max_rows=135)
        #data2= np.loadtxt(run+"He_t-1.dat", skiprows=2, max_rows=135)
        
        r0 = data_h[4,0]
        cs = data_h[-1,15]
        T0 = data_h[-1,12]
        flux0 = data_h[-1,2]*data[-1,0]**2./m0
        
        rs = G*mp/cs/cs/2.
        lambda0 = rs/r0
        #b       = 5.0e17 * (T0/300.)**0.75
        b       = 5.0e17 * (T0/1.)**0.75
        mu      = m1/m0
        phi     = flux0*kb*T0/(mu*G*mp*m0*b)
        
        phiratio  = data_he[-5,2]/data_h[-5,2]
        densratio = data_he[4,1]/data_h[4,1]
        
        return phiratio/densratio
    
#
# Returns dimensional flux and dimensionless fractionation factor (for direct comparison to simulations)
#
#
def get_fractionation(m_1,m_2, mpl, T0, r0, bfactor = 5e17):
    
    phis             = 10**np.linspace(-4,4,10000)
    #b                = bfactor * (T0/300.)**0.75
    b                = bfactor * (T0/1.)**0.75
    phi0_constant    = kb*T0/(G*mpl*mearth*m_2*amu*b) / (4.*3.141592*amu)
    redmass          = m_1*m_2/(m_1+m_2) 
    lam              = G*mpl*mearth/(kb*T0) /(r0) * amu * redmass
    xfrac_ZK86       = get_analytic_fractionation(phis, m_2/m_1, lambda0=lam, muexp=0.0) #Use isothermal Zahnle & Kasting results
    xfrac_Hu87       = get_analytic_fractionation_hunten(phis, m_2/m_1, lambda0=lam) #Use isothermal Zahnle & Kasting results

    #print("diffusion flux constant = " + repr(phi0_constant) + " dimensionless mdot=1e10: = " + repr(1e10/(12.*amu*r0**2.) * phi0_constant))
    return phis/phi0_constant, xfrac_ZK86, xfrac_Hu87



#
# Returns transit plot in gaussian mode with nonisothermal and isothermal and isothermal parker wind profiles
#
def get_transit_plots(ftitle, stitles, fsize, fdir, ssims, timesteps=[-1], cols=['black'], hefraction=[[1,1,1]], 
                      method='average', offset_x=0., offset_y=[0.],yerr=0., datapts="", T0=10000., mplanet=353., rstellar=[1.*rsolar], dotcol='black',
                     print_profile_data=0):

    w0, w1, w2, f0, f1, f2, a_ij = lines.he_3_properties()
    m_He = 4 * 1.67262192369e-27  # Helium atomic mass in kg
    wl = np.linspace(1.0827, 1.0832, 200) * 1E-6  # Wavelengths in m
    w_array = np.array([w0, w1, w2])
    f_array = np.array([f0, f1, f2])
    a_array = np.array([a_ij, a_ij, a_ij])  # This is the same for all lines in then triplet
    
    fig0 = plt.figure(figsize=fsize) 
    if ftitle != "":
        fig0.suptitle(ftitle, fontsize=20)
    ax0 = fig0.add_subplot(111)

    axes = [ax0]

    masses = np.array([5e-4, 1.,1., 4.,4.,4.,4.,4.,4.,4.,4.])*amu
    print("running sims = " + repr(ssims))
    
    for i,sim in enumerate(ssims):
        print("running sim0 = " + repr(sim))
        if 0==0 :
            print("running sim1 = " + repr(sim))
            rp = 1.
            
            dirr = fdir
            timestep = timesteps[i]
            filename = "output_"+sim+"_He11S_t" + repr(timestep)+ ".dat"    
            max_rr   = get_filelength(dirr, filename)

            data_h0 = np.loadtxt(dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_p  = np.loadtxt(dirr + "output_"+ sim +"_S1_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_e  = np.loadtxt(dirr + "output_"+ sim +"_e-_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            data_he11s  = np.loadtxt(dirr + "output_"+ sim +"_He11S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_he23s  = np.loadtxt(dirr + "output_"+ sim +"_He23S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_hep    = np.loadtxt(dirr + "output_"+ sim +"_Hep_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            data_hepp = 0.* data_hep;

            if os.path.isfile(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat"): #If O3p exist, assume that O4p also exists
                data_hepp = np.loadtxt(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            n_e  = data_e[:,1]/masses[0]; 
            n_p  = data_p[:,1]/masses[1]; 
            n_h0 = data_h0[:,1]/masses[2];

            n_he11s  = data_he11s[:,1]/masses[3];

            n_hep    = data_hep[:,1]/masses[3];
            n_hepp   = data_hepp[:,1]/masses[3];

            r_ai     = data_h0[:,0] 
            r_pl_ai  = r_ai[0] #data_h0[0,0]

            charge_ratio   = (n_p + n_hep + 2.* n_hepp)/n_e
            charge_balance = (n_p + n_hep + 2.* n_hepp - n_e)/n_e
            n_he23s        = data_he23s[:,1]/masses[3];
            n_hetot        = n_he11s + n_he23s + n_hep + n_hepp
            f_he23s        = n_he23s / n_hetot

            T_ai                 = data_h0[:,12]
            v_ai                 = data_he23s[:,11]

            #data_he23s_alt  = np.loadtxt(dirr + "output_"+ "hd209458b-helium-long" +"_He23S_t" + repr(4)+ ".dat",   skiprows=2, max_rows=max_rr)
            #v_ai_alt        = data_he23s_alt[:,11] * 1e-2

            #v_ai_alt        = data_he23s_alt[:,11] * 1e-2
            #
            # Simple, oklopcic-like calculations
            #
            
            alpha_rec         = 9.94e-11 * T_ai**(-0.668)
            alpha_rec_const   = 9.94e-11 * T0**(-0.668)
            alpha_rad         = 1.27e-4

            q_dex             = 1.02e-5 * T_ai**(-0.5) * np.exp(-9234./T_ai)
            #q_dex_const       = 2.6e-8
            q_dex_const       = 3.78e-5 * T0**(-0.69) * np.exp(-1.02e4/T0)
            n_he23s_predicted       = alpha_rec/q_dex             * n_hep  * 1e6
            n_he23s_predicted_const = alpha_rec_const/q_dex_const * n_hep  * 1e6

            # First convert everything to SI units
            
            R_pl_physical = r_pl_ai/1e2  # Planet radius in m
            r_SI          = (r_ai - r_pl_ai)/1e2 # Array of altitudes in m
            v_SI          = v_ai * 1e-2  # Velocity of the outflow in m / s
            n_he23s_SI    = n_he23s * 1E6  # Volumetric densities in 1 / m ** 3
            #planet_to_star_ratio = 0.12086
            planet_to_star_ratio = r_pl_ai / rstellar[i]

            # Set up the ray tracing. We will use a coarse 100-px grid size,
            # but we use supersampling to avoid hard pixel edges.
            impact_parameter = 0.

            flux_map3, t_depth3, r_from_planet3 = transit.draw_transit(
                planet_to_star_ratio, 
                planet_physical_radius=R_pl_physical, 
                impact_parameter=impact_parameter, 
                phase=-0.0,
                supersampling=10,
                grid_size=200)
            
            # And now we plot it just to check how the transit looks
            #plt.imshow(flux_map3, origin='lower')
            #plt.show()
            
            if hefraction[i][0] > 1e-10:
                spectrum_noniso         = transit.radiative_transfer_2d(flux_map3, r_from_planet3, r_SI, n_he23s_SI * hefraction[i][0], v_SI,  w_array, f_array, a_array, wl, T_ai, m_He, wind_broadening_method=method)
                print("Min/Max difference in spectrum = " + repr( [np.max(spectrum_noniso)- np.min(spectrum_noniso), np.max(spectrum_noniso), np.min(spectrum_noniso) ] )  )
                #ax0.plot(wl * 1E6 + offset_x, spectrum_noniso      + offset_y, ls='-',c=cols[i], label=stitles[i]+ ' x'+repr(hefraction[i][0])+' nonisothermal', lw=2.5)
                ax0.plot(wl * 1E6 + offset_x, spectrum_noniso      + offset_y[i], ls='-', c=cols[i], label=stitles[i][0], lw=2.5)
            
            if print_profile_data == 1:
                print("wl = " + repr(wl * 1E6 + offset_x))
                print("spectrum noniso = " + repr(spectrum_noniso-spectrum_noniso[0])) #APRIL29
            
            v_isotherm = analytic_wind_solution(mplanet=mplanet, mparticle=0.5, T=T0, radii=r_SI*1e2)
            v_isotherm = [v*1e-2 for v in v_isotherm]
            
            if hefraction[i][1] > 1e-10:
                spectrum_iso            = transit.radiative_transfer_2d(flux_map3, r_from_planet3, r_SI, n_he23s_SI * hefraction[i][1], v_isotherm,  w_array, f_array, a_array, wl, T0, m_He, wind_broadening_method=method)
                print("Min/Max difference in spectrum = " + repr( [ np.max(spectrum_iso)- np.min(spectrum_iso), np.max(spectrum_noniso) , np.min(spectrum_noniso) ] )  )
                #spectrum_iso            = transit.radiative_transfer_2d(flux_map3, r_from_planet3, r_SI, n_he23s_predicted_const * hefraction[i][1], v_SI,  w_array, f_array, a_array, wl, T0, m_He, wind_broadening_method=method)
                #ax0.plot(wl * 1E6 + offset_x, spectrum_iso         + offset_y, ls='--',c=cols[i], label=stitles[i]+ ' x'+repr(hefraction[i][1])+' nonisothermal T, Parker v', lw=2.5)
                ax0.plot(wl * 1E6 + offset_x, spectrum_iso         + offset_y[i], ls='--', c=cols[i], label=stitles[i][1], lw=2.5)
            
            if hefraction[i][2] > 1e-10:
                spectrum_iso2           = transit.radiative_transfer_2d(flux_map3, r_from_planet3, r_SI, n_he23s_predicted_const * hefraction[i][2], v_isotherm, w_array, f_array, a_array, wl, T0, m_He, wind_broadening_method=method)
                print("Min/Max difference in spectrum = " + repr( [ np.max(spectrum_iso2)- np.min(spectrum_iso2), np.max(spectrum_noniso), np.min(spectrum_noniso) ] ))
                ax0.plot(wl * 1E6 + offset_x, spectrum_iso2 + offset_y[i], ls=':', c=cols[i], label=stitles[i][2], lw=2.5)

            ax0.axhline(y=1.0)
            
    datadirr    = "/home/matthaus/Documents/Notebooks/aiolos/plots/data/"
    #sim     = "Salz2018_HD189733b_Heliumdata.csv"
    if datapts != "":
        data_hd = np.loadtxt(datadirr + datapts, delimiter=' ')
        #data_hd = np.loadtxt(datadirr + sim, delimiter=' ',  skiprows=2, max_rows=max_rr)
        #ax0.errorbar((data_hd[:,0] *1e-4), 1+data_hd[:,1]*1e-2, yerr=yerr, fmt="o",c='k')
        ax0.errorbar(data_hd[:,0], data_hd[:,1], yerr=yerr, fmt="o",c=dotcol)
    
    for i,spl in enumerate([ax0]):
        spl.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
        spl.tick_params(axis='x', which='major', labelsize=25, size=10,direction="in");
        spl.tick_params(axis='y', which='major', labelsize=28, size=10,direction="in");
        spl.yaxis.set_ticks_position('both')
        spl.xaxis.set_ticks_position('both')
        spl.set_xlabel('Wavelength in air ($\mu$m)', fontsize=26)
        spl.set_ylabel('Normalized flux', fontsize=26)
        #spl.set_xlim([1.0828,1.08315])
        spl.set_xlim([1.0825,1.0834])
        spl.legend(fontsize=17)

    return fig0


def plot_parameter_variation_atoms(rp, directory, sims, times, max_rrs, labels, cols, styles, title, fsize=(24.,36.), printadvectiontime=0):
    
    fig0 = plt.figure(figsize=fsize) #36,24
    #fig0.suptitle(r"$\rm 4\;m_{\oplus}, F=10^4\;erg\;cm^{-2}\;s^{-1},\;atomic\;source,\;high\;vs\;low\;hydrogen,\;cooling\;vs\;nocooling$", fontsize=20)
    #fig0.suptitle(title, fontsize=20)
    
    spl0 = fig0.add_subplot(431)
    spl1 = fig0.add_subplot(432)
    spl2 = fig0.add_subplot(433)
    spl3 = fig0.add_subplot(434)
    spl4 = fig0.add_subplot(435)
    spl5 = fig0.add_subplot(436)
    spl6 = fig0.add_subplot(437)
    spl7 = fig0.add_subplot(438)
    spl8 = fig0.add_subplot(439)

    spl9   = fig0.add_subplot(4,3,10)
    spl10  = fig0.add_subplot(4,3,11)
    spl11  = fig0.add_subplot(4,3,12)
    #spl12  = fig0.add_subplot(4,4,13)
    #spl13  = fig0.add_subplot(4,4,14)
    #spl14  = fig0.add_subplot(4,4,15)
    #spl15  = fig0.add_subplot(4,4,16)

    axes = [spl0, spl1, spl2, spl3, spl4, spl5, spl6, spl7, spl8, spl9, spl10, spl11]

    names = [
             r"$\rm Ion\,fraction\,H\;[1]$", 
             r"$\rm Ion\,fraction\,He+\;[1]$",
             r"$\rm Ion\,fraction\,He2+\;[1]$",
             r"$\rm Number\,fraction\,He2^3S/All He\;[1]$", 
             r"$\rm Number\,density\,He2^3S\;[cm^{-3}]$",
             r"$\rm Temperature\,[K]$", 
             r"$\rm Number\,density\,e^{-}(--),\;p^{+}(-)\;He^{+}(:)\;[cm^{-3}]$", #r"$\rm Number\,fraction\,He2^3S/He+\;[1]$", 
             r"$\rm Velocity\,[cm\,s^{-1}]$", #r"$\rm Velocity\,H0(-),\;He23S(--)\,[cm\,s^{-1}]$",
             r"$\rm H_{tot}\; Mass\,flux\,[m_{\oplus}\,yr^{-1}]$",
             r"$\rm Charge\; imbalance\,[n_{e-}]$", 
             r"$\rm He_{tot}\; Mass\,flux\,[m_{\oplus}\,yr^{-1}]$", 
             r"$\rm Number\,density\;n_{H,tot}(-),\,n_{He,tot}(--)\;[cm^{-3}]$",]

    masses = np.array([5e-4, 1.,1., 4.,4.,4.,4.,4.,4.,4.,4.])*amu
    
    #
    # For Documentation purposes: 10.000 outputs in this simulation correspond to 3e11s, or 10.000yrs
    #
    print(timesteps)
    for i,spl in enumerate(axes):
        #spl.set_title(names[i], fontsize=20, pad = 1)
        spl.set_ylabel(names[i], fontsize=20)
        spl.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
        spl.tick_params(axis='x', which='major', labelsize=25, size=10,direction="in");
        spl.tick_params(axis='y', which='major', labelsize=28, size=10,direction="in");
        spl.yaxis.set_ticks_position('both')
        spl.xaxis.set_ticks_position('both')
        spl.set_xlim([0.9, 10.5])

    for i,spl in enumerate([axes[-3], axes[-2], axes[-1]]):
        spl.set_xlabel(r'$\rm Radius \;[R_{\rm p}]$', fontsize=26)

    for i, spl in enumerate([spl0, spl1, spl2]): #H+,He+,He++ ion fractions
        spl.set_ylim([1e-3, 1.1])
        
    spl2.set_ylim([1e-3, 1.1]) 
        
    for i, spl in enumerate([spl3]): #H23s fraction
        spl.set_ylim([1e-10, 1e-1])
        
    #for i, spl in enumerate([spl4, spl5]): #H23s number density
    #    spl.set_ylim([1e-2, 1e+5])

    spl7.set_ylim([1e2, 1e7]) #Velocity

    for i, spl in enumerate([spl8,spl10 ]): #mass fluxes
        spl.set_ylim([1e3, 1e13])

    spl9.set_ylim([1e-1, 1e1])    
    
    spl11.set_ylim([1e+0, 1e15])   

    for i,sim in enumerate(sims):
            timestep = timesteps[i]

            max_rr =  get_filelength(dirr, "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat" ) #max_rrs[i]
            data_h0 = np.loadtxt(dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_p  = np.loadtxt(dirr + "output_"+ sim +"_S1_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_e  = np.loadtxt(dirr + "output_"+ sim +"_e-_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            data_he11s = data_h0
            data_he23s = []
            
            if os.path.isfile(dirr + "output_"+ sim +"_He11S_t" + repr(timestep)+ ".dat"):
                data_he11s = np.loadtxt(dirr + "output_"+ sim +"_He11S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            else: 
                data_he11s = 0. * data_he11s
            
            hetriplet_exists = 0
            if os.path.isfile(dirr + "output_"+ sim +"_He23S_t" + repr(timestep)+ ".dat"): #If O3p exist, assume that O4p also exists
                data_he23s  = np.loadtxt(dirr + "output_"+ sim +"_He23S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
                hetriplet_exists = 1
            else:
                data_he23s = 0 * data_he11s
            
            data_hep = data_he11s
            if os.path.isfile(dirr + "output_"+ sim +"_Hep_t" + repr(timestep)+ ".dat"): 
                data_hep = np.loadtxt(dirr + "output_"+ sim +"_Hep_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            data_hepp = 0.* data_hep;
            
            data_he21s = []
            switch_21S_exists = 0
            data_he21p = []
            switch_21P_exists = 0
            
            if os.path.isfile(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat"): #If O3p exist, assume that O4p also exists
                data_hepp = np.loadtxt(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            if os.path.isfile(dirr + "output_"+ sim +"_He21S_t" + repr(timestep)+ ".dat"): #If He21S exists
                switch_21S_exists = 1
                data_he21s = np.loadtxt(dirr + "output_"+ sim +"_He21S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            if os.path.isfile(dirr + "output_"+ sim +"_He21P_t" + repr(timestep)+ ".dat"): #If He21S exists
                switch_21P_exists = 1
                data_he21p = np.loadtxt(dirr + "output_"+ sim +"_He21P_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            data_h2 = 0.* data_h0;
            data_h2p= 0.* data_h0;
            if os.path.isfile(dirr + "output_"+ sim +"_H2_t" + repr(timestep)+ ".dat"): 
                data_hep = np.loadtxt(dirr + "output_"+ sim +"_H2_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            if os.path.isfile(dirr + "output_"+ sim +"_H2p_t" + repr(timestep)+ ".dat"): 
                data_hep = np.loadtxt(dirr + "output_"+ sim +"_H2p_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            n_e  = data_e[:,1]/masses[0]; 
            n_p  = data_p[:,1]/masses[1]; 
            n_h0 = data_h0[:,1]/masses[2];
            n_h2 = data_h2[:,1]/(2.*amu);
            n_h2p= data_h2p[:,1]/(2.*amu);

            n_he11s  = data_he11s[:,1]/masses[3];
            n_he23s  = data_he23s[:,1]/masses[3];
            n_hep    = data_hep[:,1]/masses[3];
            n_hepp   = data_hepp[:,1]/masses[3];
            
            r  = data_h0[:,0]/rp
            dr = np.diff(r)*rp

            #indexes = [i for i,val in enumerate(range(len(diag[:,0]))) if diag[i,7] < 1.]
            #if len(indexes)==0:
            #    irp = 1
            #else:
            #    irp     = indexes[0]
            #
            #r_true = r[irp]

            charge_ratio   = (n_p + n_hep + 2.* n_hepp)/n_e
            charge_balance = (n_p + n_hep + 2.* n_hepp - n_e)/n_e
            n_hetot = n_he11s + n_he23s + n_hep + n_hepp
            
            T                 = data_h0[:,12]
            Tel               = data_e[:,12]
            v                 = data_he23s[:,11]
            alpha_rec         = 9.94e-11 * T**(-0.668)
            alpha_rec_el      = 9.94e-11 * Tel**(-0.668)
            alpha_rad         = 1.27e-4
            
            q_dex             = 3.78e-5 * T**(-0.69) * np.exp(-1.015499e4/T) #1.02e-5 * T**(-0.5) * np.exp(-9234./T)
            q_dex_el          = 3.78e-5 * Tel**(-0.69) * np.exp(-1.015499e4/Tel) #1.02e-5 * Tel**(-0.5) * np.exp(-9234./Tel)
            #q_dex2            = 4.02e-6 * T**(-0.5) * np.exp(-16219./T)
            q_dex_const       = 2.6e-8
            #n_he23s_predicted = alpha_rec/alpha_rad * n_e * n_hep #Black predicts the steady state based on recombination rates and radiative decay. Actually, loss to 21P is 1e4 times faster
            n_he23s_predicted     = alpha_rec/q_dex  * n_hep
            n_he23s_predicted_tel = alpha_rec_el/q_dex_el  * n_hep
            n_he23s_predicted2    = alpha_rec/q_dex_const * n_hep
            
            colltime   = 1./(q_dex * n_e)
            advtime    = 1./(v/rp)
            recombtime = 1./(alpha_rec * n_e * n_hep/n_he23s )
            radrectime = 1./alpha_rad *n_e/n_e
             
            ionfracs    = [n_p/(n_h0 + n_p), n_hep/n_hetot, n_hepp/n_hetot, n_he23s/n_hetot, n_he23s/(n_hep), n_he23s, charge_balance]  

            #for s, spl in enumerate(axes):
            #    if s<6:
            #        spl.loglog(r, ionfracs[s], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])

            spl0.loglog(r, ionfracs[0], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])
            spl1.loglog(r, ionfracs[1], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])
            spl2.loglog(r, ionfracs[2], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])
            spl3.loglog(r, ionfracs[3], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])       
            spl4.loglog(r, ionfracs[5], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])
            
            if printadvectiontime == 1:
                spl4.loglog(r, colltime, lw=2.5, c='magenta', ls='-',  label=r"$t_{q(2^3S->2^1P)}\;[s]$")
                spl4.loglog(r, advtime, lw=2.5, c='cyan',     ls='-',  label=r"$t_{\rm adv}\;[s]$")
                spl4.loglog(r, recombtime, lw=2.5, c='orange',     ls='-',label=r"$t_{\alpha (He+,e-)}\;[s]$")
                spl4.loglog(r, radrectime, lw=2.5, c='skyblue',     ls='-',label=r"$t_{rad}\;[s]$")
                spl4.loglog(r, n_he23s_predicted2, lw=1.5, c='grey', ls=':')
                spl4.loglog(r, n_he23s_predicted, lw=1.5, c='grey', ls='--')
            elif printadvectiontime == 2:
                spl4.loglog(r, n_he23s_predicted, lw=1.25, c=cols[i], ls='--')
                spl4.loglog(r, n_he23s_predicted_tel, lw=1.25, c=cols[i], ls=':')
                spl4.loglog(r, n_he23s_predicted2, lw=1.25, c=cols[i], ls='--')
            else:
                if styles[i] != '':
                    spl4.loglog(r, n_he23s_predicted2, lw=1.5, c='grey', ls=':')
                    spl4.loglog(r, n_he23s_predicted, lw=1.5, c='grey', ls='--')
                #spl5.loglog(r, n_he23s_predicted2, lw=2.5, c='grey', ls=':', label="const q")
                #spl5.loglog(r, n_he23s_predicted, lw=2.5, c='grey', ls='--',label="variable q")
            
            spl6.loglog(r, n_e, lw=2., c=cols[i], ls='--')
            spl6.loglog(r, n_p, lw=2., c=cols[i], ls='-')
            spl6.loglog(r, n_hep, lw=2., c=cols[i], ls=':')
            spl6.loglog([0], [0], lw=2., c=cols[i], ls=styles[i], label=labels[i])

            
            spl5.loglog(r, data_he23s[:,12], lw=3., c=cols[i], ls='--')
            spl5.loglog(r, data_e[:,12],     lw=3., c=cols[i], ls=':')
            spl5.loglog(r, data_h0[:,12], lw=3., c=cols[i], ls=styles[i], label=labels[i])
            spl7.loglog(r, data_p[:,11], lw=3., c=cols[i], ls='--', label='HI')
            spl7.loglog(r, data_he23s[:,11], lw=3., c=cols[i], ls='-', label=r'$\rm 2^3S He$')
            #spl7.loglog(r, data_he11s[:,11], lw=3., c='g', ls='--', label=r'$\rm 1^1S He$')
            #spl7.loglog(r, data_hep[:,11], lw=3., c='b', ls='--', label=r'$\rm He^{+}$')

            #spl7.loglog(r, data_h0[:,15], lw=3., c='b', ls='-', label=r'$\rm c_s(HI)$')
            
            mom_h = data_h0[:,2] + data_p[:,2] + data_h2[:,2] + data_h2p[:,2]
            mom_he= data_he11s[:,2] + data_he23s[:,2] + data_hep[:,2] + data_hepp[:,2]

            spl8.loglog(r, 4.*3.1413592*r*r*rp*rp * mom_h, lw=3., c=cols[i], ls=styles[i], label=labels[i])
            spl8.loglog(r, 4.*3.1413592*r*r*rp*rp * data_e[:,2], lw=3., c=cols[i], ls='--', label=labels[i])

            #if hetriplet_exists == 1:
            spl9.loglog(r, charge_ratio,lw=3., c=cols[i], ls=styles[i], label=labels[i])
            #else:
            #    spl9.loglog(r, mom_he/mom_h,lw=3., c=cols[i], ls=styles[i], label=labels[i])
            #    spl9.set_title( r"$\rm Ratio\;of\;total\; mass\,fluxes\;\dot{m}_{He}/\dot{m}_{H}\;[1]$", fontsize=20, pad = 1)
                
            spl10.loglog(r, 4.*3.1413592*r*r*rp*rp * mom_he , lw=3., c=cols[i], ls=styles[i], label=labels[i])
            spl10.loglog(r, 4.*3.1413592*r*r*rp*rp * data_he23s[:,2] , lw=1., c=cols[i], ls=styles[i])
            
            spl11.loglog(r, n_h0+n_p , lw=3., c=cols[i], ls='-', label=labels[i])
            spl11.loglog(r, n_hetot , lw=3., c=cols[i], ls=':')

            print("ratios mdot_He/mdot_H = %1.2e" % (mom_he[-10]/mom_h[-10]) )
                            
            for s, spl in enumerate(axes):
                spl.legend(fontsize=10)
            #spl3.axvline(x=r_true, c=cols[i], ls=sstyls[i])
            #spl8.axvline(x=r_true, c=cols[i], ls=sstyls[i])
            #spl9.axvline(x=r_true, c=cols[i], ls=sstyls[i])
            
            #
            # TODO: Add from Black 1981, the simple, steady state He-triplet values n(He23S) = 7.82e-7 * T^{-0.6687} n(e-) * n(He+) cm^{-3}
            # T in Kelvin
            #

    spl4.axhline(y=1,c='k')
    #spl6.axhline(y=1,c='k')
    fig0.tight_layout()
    
    return fig0



def return_profiles(sim,timestep, max_rr=550, maskradii_percent=[]):
    
    masses = np.array([5e-4, 1.,1., 4.,4.,4.,4.,4.,4.,4.,4.])*amu

    #max_rr = 226
    #max_rr = 340
    #max_rr = 550
    dirr = "/home/matthaus/Documents/Notebooks/aiolos/plots/data/waterplanets/hdepletion/"

    #sim  = "hd189733b-helium-hdlongdeep-flattemperature-kzz1e-11_contractinginnerbound"
    #timestep = 1

    data_h0 = np.loadtxt(dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
    data_p  = np.loadtxt(dirr + "output_"+ sim +"_S1_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
    data_e  = np.loadtxt(dirr + "output_"+ sim +"_e-_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

    data_he11s  = np.loadtxt(dirr + "output_"+ sim +"_He11S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
    data_he23s  = np.loadtxt(dirr + "output_"+ sim +"_He23S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
    data_hep    = np.loadtxt(dirr + "output_"+ sim +"_Hep_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

    data_hepp = 0.* data_hep;

    if os.path.isfile(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat"): #If O3p exist, assume that O4p also exists
        data_hepp = np.loadtxt(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

    n_e  = data_e[:,1]/masses[0]; 
    n_p  = data_p[:,1]/masses[1]; 
    n_h0 = data_h0[:,1]/masses[2];

    n_he11s  = data_he11s[:,1]/masses[3];

    n_hep    = data_hep[:,1]/masses[3];
    n_hepp   = data_hepp[:,1]/masses[3];

    r_ai     = data_h0[:,0] 
    r_pl_ai  = r_ai[0] #data_h0[0,0]

    charge_ratio   = (n_p + n_hep + 2.* n_hepp)/n_e
    charge_balance = (n_p + n_hep + 2.* n_hepp - n_e)/n_e
    n_he23s              = data_he23s[:,1]/masses[3];
    n_hetot = n_he11s + n_he23s + n_hep + n_hepp
    f_he23s = n_he23s / n_hetot

    T_ai                 = data_e[:,12]
    v_ai                 = data_he23s[:,11]

    r_SI          = (r_ai - r_pl_ai)/1e2 # Array of altitudes in m
    v_SI          = v_ai * 1e-2  # Velocity of the outflow in m / s
    n_he_3_SI     = n_he23s * 1E6  # Volumetric densities in 1 / m ** 3
    
    R_pl_physical = r_pl_ai/1e2  # Planet radius in m
    
    alpha_rec         = 9.94e-11 * T_ai**(-0.668)
    alpha_rad         = 1.27e-4

    q_dex             = 1.02e-5 * T_ai**(-0.5) * np.exp(-9234./T_ai)
    q_dex2            = 4.02e-6 * T_ai**(-0.5) * np.exp(-16219./T_ai)
    q_dex_const       = 2.6e-8

    n_he23s_predicted  = alpha_rec/q_dex  * n_hep         *1e6
    n_he23s_predicted2 = alpha_rec/q_dex_const * n_hep    *1e6
    
    
    T_avg = 10000. #Tavg is the Tavg in the masked region 
    
    if len(maskradii_percent) != 0: 
        
        #Maskradii are given, mask r,n,v values and compute T_avg on the masked radius, so that T_avg can be used as isothermal value in p-winds
        
        #
        # Calculate radial cutoff from the maskradii
        #
        rmin     = np.min(r_SI); rmax = np.max(r_SI)
        rlength = len(r_SI)
        imin    = rlength*maskradii_percent[0]
        imax    = rlength*maskradii_percent[1]
        #print(" rmin = " + repr(rmin) + " rmax = "  + repr(rmax))
        #print(" len(r_SI) = " + repr(rlength))
        rangemin = 10.**( np.log10(rmin) + maskradii_percent[0] * (np.log10(rmax) - np.log10(rmin)) ) 
        rangemax = 10.**( np.log10(rmin) + maskradii_percent[1] * (np.log10(rmax) - np.log10(rmin)) ) 
        
        mask = (r_SI < rangemax) | (r_SI > rangemin)
        num_nonzeros = len(mask)
        
        num_nonzeros2 = 0
        for i,r in enumerate(r_SI):
            if i<imin or i>imax:
                n_he_3_SI[i]    = 0.
                T_ai[i]         = 0.
                v_SI[i]         = 0.
            else:
                num_nonzeros2 += 1
        
        #n_he_3_SI = np.where(mask, 0, n_he_3_SI)
        #n_he23s_predicted2
        #T_ai      = np.where(mask, 0, T_ai)
        #v_ai      = np.where(mask, 0, v_ai)
        
        T_avg = np.sum(T_ai)/float(num_nonzeros2)
        #
        #
        #
        #print("Successfully cut some data, got an average T = " + repr(T_avg))
    #else:
        #print("Nothing cut")
        
    return r_SI, n_he_3_SI, v_SI, R_pl_physical, n_he23s_predicted2, T_avg

def plot_helium_fractionation(fig, axs, rp, directory, sims, FUVs, times, max_rrs, label, col, style, markertype, markersize, plot_neutrals=0):
    
    masses = np.array([5e-4, 1.,1., 4.,4.,4.,4.,4.,4.,4.,4.])*amu
    
    n0_h         = np.asarray([])
    n0_he        = np.asarray([])
    ndot_h         = np.asarray([])
    ndot_h_neutral = np.asarray([])
    ndot_he        = np.asarray([])
    ndot_he_neutral= np.asarray([])
    charge_ratio   = np.asarray([])
    
    dirr = directory
    
    for i,sim in enumerate(sims):
            timestep = timesteps[i]

            max_rr  = get_filelength(dirr, "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat")
            data_h0 = np.loadtxt(dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_p  = np.loadtxt(dirr + "output_"+ sim +"_S1_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_e  = np.loadtxt(dirr + "output_"+ sim +"_e-_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            data_he11s  = np.loadtxt(dirr + "output_"+ sim +"_He11S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_he23s = []
            
            hetriplet_exists = 0
            if os.path.isfile(dirr + "output_"+ sim +"_He23S_t" + repr(timestep)+ ".dat"): 
                data_he23s  = np.loadtxt(dirr + "output_"+ sim +"_He23S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
                hetriplet_exists = 1
            else:
                data_he23s = 0. * data_he11s
            
            data_hep = 0.* data_he11s
            if os.path.isfile(dirr + "output_"+ sim +"_Hep_t" + repr(timestep)+ ".dat"): 
                data_hep    = np.loadtxt(dirr + "output_"+ sim +"_Hep_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            data_hepp = 0.* data_hep;
            
            data_he21s = []
            switch_21S_exists = 0
            data_he21p = []
            switch_21P_exists = 0
            
            if os.path.isfile(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat"): 
                data_hepp = np.loadtxt(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            if os.path.isfile(dirr + "output_"+ sim +"_He21S_t" + repr(timestep)+ ".dat"): 
                switch_21S_exists = 1
                data_he21s = np.loadtxt(dirr + "output_"+ sim +"_He21S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            if os.path.isfile(dirr + "output_"+ sim +"_He21P_t" + repr(timestep)+ ".dat"):
                switch_21P_exists = 1
                data_he21p = np.loadtxt(dirr + "output_"+ sim +"_He21P_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            
            n_e  = data_e[:,1]/masses[0]; 
            n_p  = data_p[:,1]/masses[1]; 
            n_h0 = data_h0[:,1]/masses[2];        

            n_he11s  = data_he11s[:,1]/masses[3]; 
            n_he23s  = data_he23s[:,1]/masses[3];
            n_hep    = data_hep[:,1]/masses[3];
            n_hepp   = data_hepp[:,1]/masses[3];
            
            r  = data_h0[:,0]/rp
            dr = np.diff(r)*rp

            nmom_h       = (data_h0[:,2] + data_p[:,2]) / masses[1]
            nmom_he      = (data_he11s[:,2] + data_hep[:,2] + data_hepp[:,2]) / masses[3]
            ndot_h_pt         = 4.*3.1413592*r*r*rp*rp * nmom_h
            ndot_e_pt         = 4.*3.1413592*r*r*rp*rp * data_e[:,2] / masses[0]
            ndot_p_pt         = 4.*3.1413592*r*r*rp*rp * data_p[:,2] / masses[1]
            ndot_h_neutral_pt = 4.*3.1413592*r*r*rp*rp * data_h0[:,2] / masses[1]
            ndot_he_pt        = 4.*3.1413592*r*r*rp*rp * nmom_he
            ndot_hep_pt       = 4.*3.1413592*r*r*rp*rp * data_hep[:,2] / masses[3]
            ndot_hepp_pt      = 4.*3.1413592*r*r*rp*rp * data_hepp[:,2] / masses[3]
            ndot_he_neutral_pt= 4.*3.1413592*r*r*rp*rp * data_he11s[:,2] / masses[3]
            
            n0_h   = np.append(n0_h,  n_h0[2]+n_p[2])
            n0_he  = np.append(n0_he, n_he11s[2]+n_hep[2] + n_hepp[2])
            ndot_h = np.append(ndot_h,  ndot_h_pt[-5])
            ndot_h_neutral = np.append(ndot_h_neutral,  ndot_h_neutral_pt[-5])
            ndot_he= np.append(ndot_he, ndot_he_pt[-5])
            ndot_he_neutral = np.append(ndot_he_neutral,  ndot_he_neutral_pt[-5])
            
            #charge_ratio_pt = (n_p + n_hep + 2.*n_hepp)/n_e
            charge_ratio_pt = (ndot_p_pt + ndot_hep_pt + 2. * ndot_hepp_pt)/ndot_e_pt
            charge_ratio    = np.append(charge_ratio, charge_ratio_pt[5])
            #print("sim =" + repr(sim) + " ndot_h = " + repr(ndot_h_pt[-5]) + " ndot_htot = " + repr(ndot_h_pt[-5]+ ndot_h_pt[-5]))
            print("ionfractions at edge: p/(p+h) =" + repr( n_p[-5]/(n_p[-5]+n_h0[-5]) ) + " hep/(allhe) = " + repr( n_hep[-5]/(n_he11s[-5] + n_hep[-5] + n_hepp[-5])  ) + " nhepp/(allhe) = " + repr(n_hepp[-5]/(n_he11s[-5] + n_hep[-5] + n_hepp[-5]) ))
    
    axs[0].loglog(FUVs, ndot_h*masses[2]         ,  lw=3. , c=col, ls=style, label=label, marker=markertype, markersize=markersize)
    
    if plot_neutrals == 1:
            #axs[0].loglog(FUVs, ndot_h_neutral*masses[2] ,  lw=1.0, c=col, ls='--', label=label+" [I]")
            #axs[1].loglog(FUVs, ndot_he_neutral*masses[3] , lw=1.0, c=col, ls='--', label=label+" [I]")
            axs[0].loglog(FUVs, ndot_h_neutral*masses[2] ,  lw=1.0, c=col, ls='--' )
            if axs.shape[0] >= 3:
                axs[1].loglog(FUVs, ndot_he_neutral*masses[3] , lw=1.0, c=col, ls='--')
            #if axs.shape[0] >= 3:
                #axs[2].loglog(FUVs, ndot_he_neutral/ndot_he,  lw=1.0, c=col, ls='--', label=label+" He [I]/[II]")
                #axs[2].loglog(FUVs, 1-ndot_h_neutral/ndot_h,  lw=1.0, c=col, ls='-', label=label+" H [II]/[H tot]")
                #axs[2].loglog(FUVs, charge_ratio,  lw=3.0, c=col, ls='-', label=label+" q^{+}/q^{-}")
    
    ###axs[1].loglog(ndot_h*masses[2], ndot_he*masses[3] , lw=3., c=col, ls=style, label=label, marker=markertype, markersize=markersize)
    ###axs[0, 1].loglog(ndot_h*masses[2], ndot_he/ndot_h, lw=3., c=col, ls=style, label=label)
    if axs.shape[0] >= 3:
        axs[1].loglog(FUVs, ndot_he*masses[3]        ,  lw=3.  , c=col, ls=style, label=label, marker=markertype, markersize=markersize)
        axs[2].loglog(FUVs, ndot_he/ndot_h/(n0_he/n0_h), lw=3., c=col, ls=style, marker=markertype, markersize=markersize, label=label)
    else:
        axs[1].loglog(FUVs, ndot_he/ndot_h/(n0_he/n0_h), lw=3., c=col, ls=style, marker=markertype, markersize=markersize, label=label)
    
    if axs.shape[0] >= 4:
        axs[3].loglog(ndot_h*masses[2], ndot_he/ndot_h / (n0_he/n0_h), lw=3., c=col, ls=style, label=label)
        #axs[3].loglog(ndot_h*masses[2], ndot_he/ndot_h / (n0_he/n0_h), lw=3., c=col, ls=style, label=label)
    

masses = [5e-4, 1.,1., 4.,4.,12.,12.,12.,16.,16.,16.]

class simobject:
    m=1.
    h=1.
    he=0.27
    c=1e-4
    o=1e-4
    cool=0.
    fuv = 1.
    kzz=1e-10
    time=-1
    name=""
    
    outnum=-1
    
    mdoth=0.
    mdothe=0.
    mdotc=0.
    mdoto=0.
    massh=0.
    massc=0.
    masso=0.
    
    def __init__(self, name, m,h,he,c,o,fuv, cool,kzz):
        self.name=name
        self.h=h
        self.c=c
        self.o=o
        self.fuv=fuv
        self.cool=cool
        self.kzz=kzz
        
    def reset_mdots(self):
        #print("resetting mdots in %dme-%1.0e-FUV%1.0e_kzz%1.0e-%dcool-%d-autogen" % (m,h,uv,kzz,cool,co))
        self.mdoth =10.**(np.random.rand(1)+10.)
        self.mdothe=10.**(np.random.rand(1)+10.)
        self.mdotc =10.**np.random.rand(1)
        self.mdoto =10.**np.random.rand(1)
        
    def read_data(self, timestep):
        dirr        = "/home/matthaus/Documents/Notebooks/aiolos/plots/data/waterplanets/hdepletion/"
        sim         = self.name
        self.outnum = timestep
        
        # open file in read mode
        count=0
        with open(dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat", 'r') as fp:
            for count, line in enumerate(fp):
                pass
        print('Total Lines', count + 1)
        max_rr = count-3
        
        print("Reading sim = " +  "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat")
        if os.path.isfile(dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat"):
            
            #max_rr = 496
            masses = [5e-4, 1.,1., 4.,4.,12.,12.,12.,16.,16.,16.]

            data_h0 = np.loadtxt(dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_p  = np.loadtxt(dirr + "output_"+ sim +"_S1_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_e  = np.loadtxt(dirr + "output_"+ sim +"_e-_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            data_he   = np.loadtxt(dirr + "output_"+ sim +"_He11S_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_hep  = np.loadtxt(dirr + "output_"+ sim +"_Hep_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            data_hepp  = np.loadtxt(dirr + "output_"+ sim +"_Hepp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            #data_c0  = np.loadtxt(dirr + "output_"+ sim +"_C0_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            #data_cp  = np.loadtxt(dirr + "output_"+ sim +"_Cp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            #data_cpp = np.loadtxt(dirr + "output_"+ sim +"_Cpp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)

            #data_o0  = np.loadtxt(dirr + "output_"+ sim +"_O0_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            #data_op  = np.loadtxt(dirr + "output_"+ sim +"_Op_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            #data_opp = np.loadtxt(dirr + "output_"+ sim +"_Opp_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            isample = max_rr-2
            
            n_e  = data_e[:,1]/masses[0]; 
            n_p  = data_p[:,1]/masses[1]; 
            n_h0 = data_h0[:,1]/masses[2];

            n_he  = data_he[:,1]/masses[3];
            n_hep = data_hep[:,1]/masses[4];

            n_c0 = 0.* data_h0[:,1]/masses[2];
            n_cp = 0.* data_h0[:,1]/masses[2];
            n_cpp = 0.* data_h0[:,1]/masses[2];

            n_o0 = 0.* data_h0[:,1]/masses[2];
            n_op = 0.* data_h0[:,1]/masses[2];
            n_opp = 0.* data_h0[:,1]/masses[2];

            r  = data_h0[isample,0]
            #dr = np.diff(r)*rp
            
            mom_h = data_h0[isample,2] + data_p[isample,2]
            mom_he= data_he[isample,2] + data_hep[isample,2] + data_hepp[isample,2]
            mom_o = 0. #data_o0[isample,2] + data_op[isample,2] + data_opp[isample,2]
            mom_c = 0. #data_c0[isample,2] + data_cp[isample,2] + data_cpp[isample,2]

            self.mdoth = 4.*3.1413592*r*r * mom_h  #* data_h0[isample,11] 
            self.mdothe= 4.*3.1413592*r*r * mom_he #* data_he[isample,11] 
            self.mdotc = 4.*3.1413592*r*r * mom_c  #* data_c0[isample,11] 
            self.mdoto = 4.*3.1413592*r*r * mom_o  #* data_o0[isample,11]
            self.massh = (data_h0[isample,20] + data_p[isample,20]) * 1e7
            self.masshe= (data_he[isample,20] + data_hep[isample,20] + data_hepp[isample,20]) * 1e7
            self.massc = 0.#(data_c0[isample,20] + data_cp[isample,20] + data_cpp[isample,20]) * 1e7 * 1e-2
            self.masso = 0.# (data_o0[isample,20] + data_op[isample,20] + data_opp[isample,20]) * 1e7 * 1e-2
        else:
            print("Not found sim = " + dirr + "output_"+ sim +"_S0_t" + repr(timestep)+ ".dat")
            self.mdoth = np.nan
            self.mdothe= np.nan
            self.mdotc = np.nan
            self.mdoto = np.nan


def compute_evolution(dataslice, times, init_masses, integrationdirection, use_diffusion_limit, Teq, rplanet, mplanet, diffprefactor=5e17):
    #print("times = " + repr(times))
    dts = np.diff(times).tolist()
    #print("dts = " + repr(dts))
    #dts.append(times[1]-times[0])
    
    #print("dts = " + repr(dts))
    #print(len(dts))
    print(len(dataslice))
    #specmasses = [init_masses]
    
    
    specmasses       = np.zeros((len(times.tolist()), len(init_masses)))
    totalmasslost    = np.zeros((len(times.tolist()), len(init_masses)))
    mdotused         = np.zeros((len(times.tolist()), len(init_masses)))
    #print(len(specmasses2[:,0]))
    #print("after init specmasses2= " + repr(specmasses2))
    specmasses[0] = np.array(init_masses)
    #print("after init specmasses2= " + repr(specmasses2))
    hdropcounter = 0
    cnt          = 0
    incr         = +1
    stepcnt      = 0
    switchtimes  = []
    if integrationdirection == -1:
        cnt  = dataslice.size
        incr = -1
    
    slicesize = dataslice.size
    mdotipp = np.zeros((slicesize, len(init_masses)+1)) #matrix entry: (mh, mdoth, mdotc, mdoto) times mh timesteps
    for j,elm in enumerate(dataslice):
        mdotipp[slicesize - j-1,0] = init_masses[0] * 10.**(-j)
        mdotipp[slicesize - j-1,1] = elm.mdoth
        mdotipp[slicesize - j-1,2] = elm.mdothe
        mdotipp[slicesize - j-1,3] = elm.mdoto
    
    mdot_diffs = [] #array to save diffusion limited values at key points in the evolution
    mdot_diff  = -1 #tmp value of diffusion limited flux
    first_mdotdiffswitch = 0
    
    for t,time in enumerate(dts):

        diffusion_limited_switch = 0
        
        #
        # Step 1: interpolate mass loss rates on simulation data grid
        #
        if cnt < len(dataslice) - 1:
            mdotused[t+1,0] = 10.**np.interp(np.log10(specmasses[t,0]), np.log10(mdotipp[:,0]) , np.log10(mdotipp[:,1]))
            mdotused[t+1,1] = 10.**np.interp(np.log10(specmasses[t,0]), np.log10(mdotipp[:,0]) , np.log10(mdotipp[:,2]))   #np.interp(specmasses2[t,0], mdotipp[:,0] , mdotipp[:,3])#dataslice[cnt].mdotc
            mdotused[t+1,2] = 10.**np.interp(np.log10(specmasses[t,0]), np.log10(mdotipp[:,0]) , np.log10(mdotipp[:,3]))
            # 
            # Typical intpolatoin with numpy: #numpy.interp(x, xp, yp), where x is the free value between some of the points xp, yp
            ## Explanation for this interpolation: as the hydrogen mass controls the timesteps, the x-axis is the hydrogen mass, and the function yp is each species own escape rate
            ## i.e. we are interpolating something like  numpy.interp(m_h_current, m_h_points, mdot_whicheverspecies_points)
            #
        else:
            mdotused[t+1,0] = dataslice[cnt].mdoth
            mdotused[t+1,1] = dataslice[cnt].mdotc
            mdotused[t+1,2] = dataslice[cnt].mdoto
        
        #
        # Step 2: Check if the hydrogen mass is lower than the mass in metals - if so and if the diffusion limited switch is set, then replace the measured hydrogen escape rate with the diffusion limited one
        #
        
        if( (specmasses[t][0] < specmasses[t][1]/12. or  specmasses[t][0] < specmasses[t][2]/16.) and use_diffusion_limit == 1):
            #b           = 2.7e17*(Teq)**0.75 #From Catling&Kasting book Chapt 5, p. 145
            b           = diffprefactor*(Teq)**0.75 #
            mu          = (specmasses[t,0] + specmasses[t,1] + specmasses[t,2] ) / (specmasses[t,0]/1. + specmasses[t,1]/12. + specmasses[t,2]/16.) * amu
            n_div_ntot  = specmasses[t,0]/1./(specmasses[t,0]/1. + specmasses[t,1]/12. + specmasses[t,2]/16.)
            H           = kb*Teq/mu/(G*mplanet/rplanet**2)
            diffnumflux = b/H * n_div_ntot
            mdot_diff = 4.*3.141592*rplanet**2 * amu * diffnumflux
            
            diffusion_limited_switch = 1
            mdotused[t+1,0]          = mdot_diff
            #mdotused[t+1,1] = 10.**np.interp(np.log10(specmasses[t,0]), np.log10(mdotipp[:,0]) , np.log10(mdotipp[:,2]))   #np.interp(specmasses2[t,0], mdotipp[:,0] , mdotipp[:,3])#dataslice[cnt].mdotc
            #mdotused[t+1,2] = 10.**np.interp(np.log10(specmasses[t,0]), np.log10(mdotipp[:,0]) , np.log10(mdotipp[:,3]))
            
            if first_mdotdiffswitch == 0:
                print("diffusion limit parts at first diffusion: H[km] = " + repr(H/1e5) + " b = " + repr(b) + " numdensratio = " +  repr(n_div_ntot) + " log10(mdot) = " + repr(np.log10(mdot_diff)))
                first_mdotdiffswitch = 1

            
            if(tempmasses[0] < init_masses[0] * 10**(hdropcounter-1.)):
                print("diffusion limit parts: H[km] = " + repr(H/1e5) + " b = " + repr(b) + " numdensratio = " +  repr(n_div_ntot))
            
        
        #
        # Step 3: Compute the new masses with the mass-loss rates, either implicityl or explicitely
        #
        tempmasses     = specmasses[t]
        #cc             = dataslice[cnt].mdoth / specmasses[t,0] 
        if 0==0:
            tempmasses[0] = specmasses[t,0] / (1. + dts[t] * mdotused[t+1,0] / specmasses[t,0]  )# tempmasses[0] = specmasses[t,0] / (1. + dts[t] * dataslice[cnt].mdoth / specmasses[t,0]  ) #specmasses2[t,0] / (1. + dts[t] * cc)
            tempmasses[1] = specmasses[t,1] / (1. + dts[t] * mdotused[t+1,1] / specmasses[t,1]  ) 
            tempmasses[2] = specmasses[t,2] / (1. + dts[t] * mdotused[t+1,2] / specmasses[t,2]  ) 
        else:
            tempmasses[0] -= dts[t] * mdotused[t+1,0] #specmasses2[t,0] / (1. + dts[t] * cc)
            tempmasses[1] -= dts[t] * mdotused[t+1,1]
            tempmasses[2] -= dts[t] * mdotused[t+1,2]
        
        #print("time index = " + repr(t) + ", current masses = " + repr(tempmasses) + " mdot_h = " + repr(dataslice[cnt].mdoth))
        
        #
        # Step 4: Use the time to determin when to switch on data, with dataslice[cnt].time < time as comparison to control the escape rates switch
        #
        #
        #if(tempmasses[0] < init_masses[0] * 10**(hdropcounter-1.)):
        if(time > dataslice[cnt].time):
            
            if(cnt < len(dataslice)-1):
                cnt          += incr
                hdropcounter -= 1
                switchtimes.append(times[t])
                print("switching mass loss rates at time index = " + repr(t) + ", current masses = " + repr(tempmasses) + " time/Myr and guiding time/Myr = " + repr([time/3e13, dataslice[cnt].time/3e13]))
                
                if diffusion_limited_switch == 1:
                    mdot_diffs.append([times[t], mdot_diff])
        #else:
            #print("time index = " + repr(t) + ", current masses = " + repr(tempmasses))
        
        #specmasses = np.append([specmasses],[tempmasses], axis=0)
        #specmasses += [tempmasses]
        specmasses[t+1] = np.array(tempmasses)
        
        if 0==1:
            totalmasslost[t+1,0] = dts[t] * mdotused[t+1,0]
            totalmasslost[t+1,1] = dts[t] * mdotused[t+1,1]
            totalmasslost[t+1,2] = dts[t] * mdotused[t+1,2]
        else:
            totalmasslost[t+1,0] = init_masses[0] - specmasses[t+1,0]#totalmasslost[t,0] + dts[t] * dataslice[cnt].mdoth
            totalmasslost[t+1,1] = init_masses[1] - specmasses[t+1,1] #dts[t] * dataslice[cnt].mdotc
            totalmasslost[t+1,2] = init_masses[2] - specmasses[t+1,2] #totalmasslost[t,2] + dts[t] * dataslice[cnt].mdoto
        
        stepcnt += 1
        if stepcnt > 10:
            break
        #for s,spc in enumerate(init_masses):
        #    totalmasslost[t+1,s] = (specmasses2[t+1,s]-specmasses2[t,s]) #totalmasslost[t,s] + 
        
        print("time index = " + repr(t) + ", current specmasses = " + repr(specmasses[t+1]) + " mdots = " + repr(mdotused[t+1,:]))
        #print("time index = " + repr(t) + ", current tempmasses = " + repr(tempmasses))

        
    #specmasses = np.array(specmasses)
    #print("tempmasses = " +  repr(np.array(specmasses)[:,0] ))
    #print("specmasses2 = " +  repr(specmasses2 ))
    
    return specmasses, totalmasslost, mdotused, switchtimes, np.asarray(mdot_diffs)

#
# Function to generate two plot matrices, and separate single plot cutouts "passportphotos" to put into a/the paper
#
#     Matrix1.size = 4 x n_iontypes, every iontype is something for which ionfraction or neutralfraction and total mass flux is well-defined
#                    Rows in order: [ion or neutral fraction], [mixing ratio, total or split w/ components], [mass fluxes split or w/ components], [specials, like C/O ratio, separate component velocities]
#     Matrix2.size = 2 x 4, average quantities for the simulation: [T, v_mean (ions and neutrals), P, n ], [H, mu, -Phi, charge]
#     Passportphotos: two lists with 1 or 0 for each element of Matrix1 and Matrix2
#
# Arguments:
#                      rp: planet radius #superfluous? should be self-computed from radius array
#               directory:
#    mixingratio_switches:
#                 species:  list of length n_species, containing [name, mass, charge]
#         iontypes/groups:  list of [species name, indexdummy], as in the simulation that belong together. For every iontype j we generate a column in Matrix1 and the neutralfraction is computed as 
#                           iontypes[j][0] / sum_k iontypes[j][k] or the ion fraction as 1 - that number.
#                           indexdummy is determined in the plot function, to avoid errors when switching up simulation setups
#                           TODO: How to check for higher ionic states?
#
# ... missing documentation ...
#   plot_group_massloss_ratio: Two numbers, a and b. In the top right plot of the bottom figure, we will plot the ratio of the total mass-loss rates of the groups a and group b.
#                              e.g. if there are five groups, group[2] is the total H continaing atoms, molecules and ions and group[4] is the total Helium, and one wants the
#                              ratio of total Helium to total Hydrogen, to determine the fractionation regime, then pass plot_group_massloss_ratio = [4,2]
def plot_parameter_variation_kzzMultiatoms(directory, sims, times, labels, cols, styles, title,  
                                           species, iontypes, ionnames_total, pressure_on_xaxis, ionfrac_switches, mixingratio_switches, massflux_switches, passportfotos,
                                           fsize=(24.,36.), fsize2=(24.,36.), rplanet=-1, plot_group_massloss_ratio = [-1,0]):
    
    num_ions = len(iontypes)
    if(num_ions < 1):
        print("ERROR: num_ions has to be >= 1!   num_ions = " + repr(num_ions))
    #else:
    #    print(" num_ions = " + repr(num_ions))
    
    fig0, axs0 = plt.subplots(4, num_ions,figsize=fsize)
    fig1, axs1 = plt.subplots(2,4,figsize=fsize2)
    
    names0 = [
             r"$\rm Molecular\;fraction,H_2/H(--),\; Ion\,fraction\,H(-)\;[1]$", 
             r"$\rm Ion\,fraction\,He\;[1]$",
             r"$\rm Ion\,fraction\,C\;[1]$",
             r"$\rm Ion\,fraction\,O\;[1]$",
             r"$\rm H_2^{+}(--),H_3^{+}(-)\;Number\;densities\;[cm^{-3}]$",
             r"$\rm He2^3S\; Number\;density\;[cm^{-3}]$", 
             r"$\rm C/O\,Number\,ratio\;[1]$",
             r"$\rm Charge\; balance\,[n_{e-}]$", 
             r"$\rm H_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",
             r"$\rm He_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",
             r"$\rm C_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",
             r"$\rm O_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",]
    
    names1 = [r"$\rm Temperature\,[K]$", 
             r"$\rm \bar{v}_{neutrals}, \bar{v}_{ions}\,[cm\,s^{-1}]$",
             #r"$\rm Pressures:\; \;P_{e-}(--),\, P_{tot}(-)\,[bar]$",
             r"$\rm Pressures:\; \; P_{tot}(-),\;P_{He\,tot}(:),\,[bar]$",
             #r"$\rm Total\; number\; density\; [cm^{-3}]$", 
             r"$\rm Mass\, loss\;ratios\;\dot{M}_{He,\,tot}/\dot{M}_{H,\,tot}$", 
             r"$\rm Scale\;Height \;[R_p]$",
             r"$\rm Mean-molecular\;weight \;[amu]$",
             r"$\rm Grav.\,acceleration\;[ cm/s^2 ]$",
             r"$\rm Charge\,balance\;[ n_{e^-}/\Sigma_i q_{i} n_{i} ]$", ]
        
    #New rows in 4x4 plot matrix
    #Ion and molecular state H, He, C, O
    #mixing ratio of a representative plot? Temperature, velocity 
    #Charge balance, C/O ratio, Helium triplet, H3+
    #mdots H, He, C O

    #masses = np.array([5e-4, 1.,1., 4.,4.,4.,4.,4.,4.,4.,4.])*amu
    
    print(timesteps)
    for i,ax in enumerate(axs1.flat):
        #ax.set_title(names[i], fontsize=20, pad = 1)
        ax.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
        ax.tick_params(axis='x', which='major', labelsize=20, size=10,direction="in");
        ax.tick_params(axis='y', which='major', labelsize=20, size=10,direction="in");
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
    for i,ax in enumerate(axs0.flat):
        #ax.set_title(names[i], fontsize=20, pad = 1)
        ax.tick_params(axis='both', which='minor', labelsize=14, size=7,direction="in")
        ax.tick_params(axis='x', which='major', labelsize=16, size=10,direction="in");
        ax.tick_params(axis='y', which='major', labelsize=16, size=10,direction="in");
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    
    #for i,ax in enumerate(axs0.flat):
    #    ax.set_title(names0[i], fontsize=20, pad = 1)
    for i,ax in enumerate(axs1.flat):
        #ax.set_title(names1[i], fontsize=20, pad = 1)
        ax.set_ylabel(names1[i], fontsize=18)
    
    for t,iontype in enumerate(iontypes):
                axs0[3,t].set_ylim([1e0,1e6])
    
    masslosses_H = []
    masslosses_He =[]
    
    for i,sim in enumerate(sims):
            timestep = timesteps[i]
            
            timestep_mo_exists = 0
            timestep_po_exists = 0
            
            max_rr  = get_filelength(dirr, "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat" ) #max_rrs[i]
            data_h0 = np.loadtxt(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            if rplanet < 0:
                rp = data_h0[2,0]
            else:
                rp = rplanet
            
            if os.path.isfile(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep-1)+ ".dat"):
                timestep_mo_exists = 1
            if os.path.isfile(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep+1)+ ".dat"):
                timestep_po_exists = 1
            
            defaultr    = 0. * data_h0[:,0]        #An array of 0 elements with the correct radial extent, hence detfault"r"
            defaultdata = 0. * np.asarray(data_h0)
            
            #all_species_strings = ["e-"]   + h_species_strings + he_species_strings + c_species_strings + o_species_strings
            #all_species_datas   = [data_e] + h_species_datas   + he_species_datas   + c_species_datas   + o_species_datas
            #all_species_charges = [-1]     + h_charges         + he_charges         +  c_charges        + o_charges
            
            
            #
            # Setup global array of all species data
            #
            all_species_datas = []
            for s, spc in enumerate(species):
                file = "output_"+ sim +"_"+species[s][0]+"_t" + repr(timestep)+ ".dat"
                #print( "loading speciesdata = " + repr( file ))
                
                tmpdata = []
                if os.path.isfile(dirr + file):
                    tmpdata               = np.loadtxt(dirr + file,   skiprows=2, max_rows=max_rr)
                else:
                    tmpdata = 0. * np.loadtxt(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
                
                #all_species_datas = np.append(all_species_datas, tmpdata)
                if s == 0:
                    all_species_datas = [tmpdata]
                else:
                    #all_species_datas = np.append(all_species_datas, tmpdata)
                    all_species_datas = all_species_datas + [tmpdata]
            
            #
            # Set x-axis variable
            #
            r  = all_species_datas[0][:,0]/rp
            if pressure_on_xaxis == 1:
                r  = total_p/1e6
            #dr = np.diff(r)*rp
            
            total_n               = defaultr
            total_p               = defaultr
            total_dens_neutrals   = defaultr
            total_dens_ions       = defaultr
            total_charge          = defaultr
            avg_v_ions            = defaultr
            avg_v_neutrals        = defaultr
            avg_T                 = defaultr
            
            #Helium gets special treatment in this notebook
            total_p_heliums       = defaultr
            mmw                   = defaultr
            scaleH                = defaultr
            grav                  = defaultr
            
            print("Done with global array setup")
            #print("all_species_datas = " + repr(all_species_datas))
            #print("all_species_datas[0] = " + repr(all_species_datas[0]))
            print(" len(allspecies)" + repr(len(all_species_datas)))
            #print("all_species_datas[0][0,0] = " + repr(all_species_datas[0][0,0]))
            
            #
            # Compute total and averaged quantities over all species
            #
            for s, spc in enumerate(all_species_datas):
                total_n      = np.add(total_n, spc[:,4])
                total_charge = np.add(total_charge, species[s][2] * spc[:,4]) 
                total_p      = np.add(total_p, spc[:,10])
                #mmw          = np.add(mmw, species[s][1] * spc[:,4])
                mmw          = np.add(mmw, spc[:,1]/amu)
                avg_T        = np.add(avg_T, spc[:,4] * spc[:,12])
                
                if( species[s][2]== 0 ):
                    avg_v_neutrals      = np.add(avg_v_neutrals, spc[:,2])    
                    total_dens_neutrals = np.add(total_dens_neutrals, spc[:,1])
                else:
                    avg_v_ions     = np.add(avg_v_ions, spc[:,2])
                    total_dens_ions= np.add(total_dens_ions, spc[:,1])
                    
                if np.abs(1. - species[s][1]/4.) < 0.1: #Found helium species
                    total_p_heliums = np.add(total_p_heliums, spc[:,10])
                
            charge_ratio    = 1. + np.divide(total_charge , total_n)
            avg_v_neutrals /= total_dens_neutrals
            avg_v_ions     /= total_dens_ions
            
            avg_T          /= total_n
            mmw            /= total_n
            
            grav   = -all_species_datas[0][:,19] / (all_species_datas[0][:,0])
            scaleH = kb * avg_T / (mmw * amu * grav)
            
            #
            # For all iontypes and all ions per iontype, find the index in the all_species_data array
            #
            # e.g. iontypes with two types, for H and He would look like
            # iontypes = [     [["H0", 0] , ["p+", 0 ]],
            #                  [["He", 0] , ["He+",0] ],  
            #            ]
            for t, iontype in enumerate(iontypes): 
                for s,spc in enumerate(iontype):    #Every iontype is a list of [name, 0], where 0 needs to be replaced by the correct index in the species array
                    for ss,sspc in enumerate(species):  #Every sspc in species is a list [name, mass, charge]
                        if sspc[0] == spc[0]:
                            #index found
                            spc[1] = ss
            
            #print("After the index finding routine, iontypes = " + repr(iontypes))
            #
            # Compute total quantities for all iontypes
            #
            timestep_mo_exists = 0
            iontype_data = [] # Will be set of [ion_fraction, mixing_ratio, mass_fluxes] per iontype
            for t, iontype in enumerate(iontypes):
                t_dens     = defaultr
                t_massdens = defaultr
                
                t_massflux = defaultr
                t_massflux_po = defaultr
                t_massflux_mo = defaultr
                
                for s,spc in enumerate(iontype):
                    t_dens     = np.add(t_dens,     all_species_datas[ spc[1] ][:,4]  )
                    t_massdens = np.add(t_massdens, all_species_datas[ spc[1] ][:,1]  )
                    t_massflux = np.add(t_massflux, all_species_datas[ spc[1] ][:,5]  )
                    
                    if(timestep_po_exists):
                        tmp_data        = np.loadtxt(dirr + "output_"+ sim +"_"+spc[0]+"_t" + repr(timestep+1)+ ".dat",   skiprows=2, max_rows=max_rr)
                        t_massflux_po   = np.add(t_massflux_po, tmp_data[:,5])
                    if(timestep_mo_exists):
                        tmp_data        = np.loadtxt(dirr + "output_"+ sim +"_"+spc[0]+"_t" + repr(timestep-1)+ ".dat",   skiprows=2, max_rows=max_rr)
                        t_massflux_mo   = np.add(t_massflux_mo, tmp_data[:,5])
                
                t_ionfrac         = all_species_datas[ iontype[0][1] ][:,4] / t_dens
                t_mixingratio     = all_species_datas[ iontype[0][1] ][:,4] / total_n
                t_mixingratio_tot = t_dens / total_n
                
                if t == 0:
                    iontype_data =                 [[t_ionfrac, t_mixingratio, t_massflux, t_massflux_mo, t_massflux_po, t_mixingratio_tot, t_massdens] ]
                else:
                    iontype_data = iontype_data +  [[t_ionfrac, t_mixingratio, t_massflux, t_massflux_mo, t_massflux_po, t_mixingratio_tot, t_massdens]]
            
            print("len(iontype data) =" +  repr(len(iontype_data)))
            #print("iontype+data = " + repr(iontype_data))
            #FEB18th
            #
            # Start plotting
            #
            axs0[0,0].set_ylabel(r"$\rm Species \; ratio\;\;  n_s/\sum_{group} n_s$",fontsize=18)
            axs0[1,0].set_ylabel(r"$\rm Mixing\; ratio\;\;  n_s/\sum_{all} n_s$",fontsize=18)
            axs0[2,0].set_ylabel("Mass-loss rate [g/s]",fontsize=18)
            axs0[3,0].set_ylabel("Velocities [cm/s]",fontsize=18)
            
            fracs = []
            
            groupstyles = ['-','--',':','-.']
            for t,iontype in enumerate(iontypes):
                axs0[0,t].set_title(ionnames_total[t],fontsize=18)
                axs0[0,t].loglog(r, iontype_data[t][0], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])      #species fraction among total ion
                axs0[1,t].loglog(r, iontype_data[t][1], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])      #species mixing ratio
                axs0[1,t].loglog(r, iontype_data[t][5], lw=1., c=cols[i], ls=':')            #total mixing ratio
                axs0[2,t].loglog(r, iontype_data[t][2], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])      #mass flux
                
                if(timestep_mo_exists):
                    axs0[2,t].fill_between(r, iontype_data[t][2], iontype_data[t][3], color=cols[i], alpha=0.2, hatch='-') 
                if(timestep_po_exists):
                    axs0[2,t].fill_between(r, iontype_data[t][2], iontype_data[t][4], color=cols[i], alpha=0.2, hatch='-') 
                
                fluxratio_at_bound =  iontype_data[t][2][-1] /  iontype_data[0][2][-1]
                densratio_at_bound =  iontype_data[t][6][2]  /  iontype_data[0][6][2]
                
                if t==2:
                    masslosses_H = np.append(masslosses_H,   iontype_data[t][2][-1])
                if t==3:
                    masslosses_He = np.append(masslosses_He, iontype_data[t][2][-1])
                
                print("H fractionation factor: " + repr(fluxratio_at_bound/densratio_at_bound) + " fluxratio, densratio = " + repr([fluxratio_at_bound, densratio_at_bound]) )
                fracs = np.append(fracs, fluxratio_at_bound/densratio_at_bound)
                
                for s,spc in enumerate(iontype):
                    #axs0[3,t].loglog(r,     all_species_datas[ spc[1] ][:, 11], c=cols[i], ls=groupstyles[s], lw=2., label= spc if i==0 else ""  )
                    axs0[3,t].loglog(r,     all_species_datas[ spc[1] ][:, 11], c=cols[i], ls=groupstyles[s], lw=2., label= species[spc[1]][3] if i==0 else ""  )
                
            print("All fractionation factors: " + repr(fracs))
            #print("All H ")
            
            #axs1[0,0].loglog(r, all_species_datas[0][:,12], lw=2., c=cols[i], ls=sstyls[i], label=labels[i])      #C/O ratio in flow
            axs1[0,0].loglog(r, avg_T, lw=2., c=cols[i], ls=sstyls[i], label=labels[i])      #C/O ratio in flow
            axs1[0,1].loglog(r, avg_v_neutrals, lw=2., c=cols[i], ls=sstyls[i], label=labels[i])
            axs1[0,1].loglog(r, avg_v_ions, lw=1., c=cols[i], ls=sstyls[i], label='')
            axs1[0,2].loglog(r, total_p/1e6,         lw=3., c=cols[i], ls=sstyls[i], label=labels[i])     
            axs1[0,2].loglog(r, total_p_heliums/1e6,         lw=1., c=cols[i], ls=':')     
            #axs1[0,3].loglog(r, total_n,             lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            
            axs1[0,3].loglog(r, iontype_data[plot_group_massloss_ratio[0]][2]/iontype_data[plot_group_massloss_ratio[1]][2],             lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            #axs1[0,3].loglog(r, iontype_data[plot_group_massloss_ratio[0]][2],             lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            #axs1[0,3].loglog(r, iontype_data[plot_group_massloss_ratio[1]][2],             lw=1, c=cols[i], ls='--', label=labels[i])    
            
            #axs1[1,2].loglog(r, data_e[:,10]/1e6, lw=1., c=cols[i], ls='--', label='')
            axs1[1,0].loglog(r, scaleH / rp, lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            axs1[1,1].semilogx(r, mmw, lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            axs1[1,2].loglog(r, grav, lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            axs1[1,3].loglog(r, charge_ratio, lw=1.5, c=cols[i], ls=sstyls[i], label=labels[i])
            
            
            for a, ax in enumerate(axs0.flat):
                ax.legend(fontsize=11)
            for a, ax in enumerate(axs0[3,:]):
                ax.legend(fontsize=14)
            for a, ax in enumerate(axs1.flat):
                ax.legend(fontsize=11)
    
    print(" All mdot H = " + repr(masslosses_H))
    print(" All mdot He = " + repr(masslosses_He))
    
    xlabelstring = r"$\rm Radius \; [R_p]$"
    if rplanet > 0: #Assuming the radius given is in earth radii
        xlabelstring = r"$\rm Radius \; [R_{\oplus}]$"
        axs1[1,0].set_ylabel(r"$\rm Scale\;Height \; [R_{\oplus}]$", fontsize=16)
    if pressure_on_xaxis == 1:
        xlabelstring = r"$\rm Pressure \; [bar]$"
    
    for a, ax in enumerate(axs0[3,:]):
                ax.set_xlabel(xlabelstring,fontsize=16)
    for a, ax in enumerate(axs1[1,:]):
                ax.set_xlabel(xlabelstring,fontsize=16)
            
    #axs[1,2].axhline(y=1,c='k')
    #spl5.axhline(y=1,c='k')
    #FEB4
    
    fig0.tight_layout()
    fig1.tight_layout()
    
    return fig0, fig1
 
#
# Function to generate two plot matrices, and separate single plot cutouts "passportphotos" to put into a/the paper
#
#     Matrix1.size = 4 x n_iontypes, every iontype is something for which ionfraction or neutralfraction and total mass flux is well-defined
#                    Rows in order: [ion or neutral fraction], [mixing ratio, total or split w/ components], [mass fluxes split or w/ components], [specials, like C/O ratio, separate component velocities]
#     Matrix2.size = 2 x 4, average quantities for the simulation: [T, v_mean (ions and neutrals), P, n ], [H, mu, -Phi, charge]
#     Passportphotos: two lists with 1 or 0 for each element of Matrix1 and Matrix2
#
# Arguments:
#                      rp: planet radius #superfluous? should be self-computed from radius array
#               directory:
#    mixingratio_switches:
#                 species:  list of length n_species, containing [name, mass, charge]
#         iontypes/groups:  list of [species name, indexdummy], as in the simulation that belong together. For every iontype j we generate a column in Matrix1 and the neutralfraction is computed as 
#                           iontypes[j][0] / sum_k iontypes[j][k] or the ion fraction as 1 - that number.
#                           indexdummy is determined in the plot function, to avoid errors when switching up simulation setups
#                           TODO: How to check for higher ionic states?
#
# ... missing documentation ...
#   plot_group_massloss_ratio: Two numbers, a and b. In the top right plot of the bottom figure, we will plot the ratio of the total mass-loss rates of the groups a and group b.
#                              e.g. if there are five groups, group[2] is the total H continaing atoms, molecules and ions and group[4] is the total Helium, and one wants the
#                              ratio of total Helium to total Hydrogen, to determine the fractionation regime, then pass plot_group_massloss_ratio = [4,2]
def plot_parameter_variation_kzzMultiatoms2(directory, sims, times, labels, cols, styles, title,  
                                           species, iontypes, ionnames_total, pressure_on_xaxis, ionfrac_switches, mixingratio_switches, massflux_switches, passportfotos,
                                           fsize=(24.,36.), fsize2=(24.,36.), fsize3=(16,11), rplanet=-1, plot_group_massloss_ratio = [-1,0], snapshot=[0,0], pradii=[14], 
                                            outlabel='default', plot_isobars_mask=[[],[]], plot_isobars_labels=[] ):
    
    num_ions = len(iontypes)
    if(num_ions < 1):
        print("ERROR: num_ions has to be >= 1!   num_ions = " + repr(num_ions))
    #else:
    #    print(" num_ions = " + repr(num_ions))
    
    fig0, axs0 = plt.subplots(4, num_ions,figsize=fsize)
    fig1, axs1 = plt.subplots(2,4,figsize=fsize2)
    fig2, axs2 = plt.subplots(len(snapshot),2,figsize=fsize3)
    fig3, axs3 = plt.subplots(1,2,figsize=(20,7))
    ax77 = axs1[1,3].twiny()
    
    names0 = [
             r"$\rm Molecular\;fraction,H_2/H(--),\; Ion\,fraction\,H(-)\;[1]$", 
             r"$\rm Ion\,fraction\,He\;[1]$",
             r"$\rm Ion\,fraction\,C\;[1]$",
             r"$\rm Ion\,fraction\,O\;[1]$",
             r"$\rm H_2^{+}(--),H_3^{+}(-)\;Number\;densities\;[cm^{-3}]$",
             r"$\rm He2^3S\; Number\;density\;[cm^{-3}]$", 
             r"$\rm C/O\,Number\,ratio\;[1]$",
             r"$\rm Charge\; balance\,[n_{e-}]$", 
             r"$\rm H_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",
             r"$\rm He_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",
             r"$\rm C_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",
             r"$\rm O_{tot}(-)\; Mass\,flux\,[g\;s^{-1}]$",]
    
    names1 = [r"$\rm Temperature\,[K]$", 
             r"$\rm \bar{v}_{neutrals}, \bar{v}_{ions}\,[cm\,s^{-1}]$",
             #r"$\rm Pressures:\; \;P_{e-}(--),\, P_{tot}(-)\,[bar]$",
             r"$\rm Pressures:\; \; P_{tot}(-),\;P_{He\,tot}(:),\,[bar]$",
             #r"$\rm Total\; number\; density\; [cm^{-3}]$", 
             r"$\rm Fractionation\, factor\;\dot{M}_{He,\,tot}/\dot{M}_{H,\,tot}\times \left( \rho_{He}/\rho_{H,\,tot}(r=r_p) \right)^{-1}$", 
             r"$\rm Scale\;Height \;[R_p]$",
             r"$\rm Mean-molecular\;weight \;[amu]$",
             r"$\rm Grav.\,acceleration\;[ cm/s^2 ]$",
             r"$\rm Total\;Pressure\;\,[bar]$",
              #r"$\rm Charge\,balance\;[ n_{e^-}/\Sigma_i q_{i} n_{i} ]$", 
             ]
    
    #New rows in 4x4 plot matrix
    #Ion and molecular state H, He, C, O
    #mixing ratio of a representative plot? Temperature, velocity 
    #Charge balance, C/O ratio, Helium triplet, H3+
    #mdots H, He, C O

    #masses = np.array([5e-4, 1.,1., 4.,4.,4.,4.,4.,4.,4.,4.])*amu
    
    print(timesteps)
    for i,ax in enumerate(axs0.flat):
        #ax.set_title(names[i], fontsize=20, pad = 1)
        ax.tick_params(axis='both', which='minor', labelsize=14, size=7,direction="in")
        ax.tick_params(axis='x', which='major', labelsize=16, size=10,direction="in");
        ax.tick_params(axis='y', which='major', labelsize=16, size=10,direction="in");
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    
    for iterlist in [axs1.flat, [ax77], axs2.flat, axs3.flat]:
        for i,ax in enumerate(iterlist):
            #ax.set_title(names[i], fontsize=20, pad = 1)
            ax.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
            ax.tick_params(axis='x', which='major', labelsize=20, size=10,direction="in");
            ax.tick_params(axis='y', which='major', labelsize=20, size=10,direction="in");
            ax.yaxis.set_ticks_position('both')
            ax.xaxis.set_ticks_position('both')
    
    for i,ax in enumerate([ax77]):
        ax.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
        ax.tick_params(axis='x', which='major', labelsize=20, size=10,direction="in");
        ax.tick_params(axis='y', which='major', labelsize=20, size=10,direction="in");
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
    for i,ax in enumerate(axs2.flat):
        ax.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
        ax.tick_params(axis='x', which='major', labelsize=20, size=10,direction="in");
        ax.tick_params(axis='y', which='major', labelsize=20, size=10,direction="in");
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        
    for i,ax in enumerate(axs3.flat):
        ax.tick_params(axis='both', which='minor', labelsize=18, size=7,direction="in")
        ax.tick_params(axis='x', which='major', labelsize=20, size=10,direction="in");
        ax.tick_params(axis='y', which='major', labelsize=20, size=10,direction="in");
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    
    for i,ax in enumerate(axs1.flat):
        ax.set_ylabel(names1[i], fontsize=18)
        
    for i,ax in enumerate(axs2.flat):
        ax.set_ylabel("Mixing ratio", fontsize=18)
    
    for t,iontype in enumerate(iontypes):
                axs0[3,t].set_ylim([1e0,1e7])
    
    ax77.set_ylim([1e3,1e-15])
    ax77.set_xlim([1e-2,1e7])
    axs1[1,3].set_xlim([450,20000.])
    axs1[1,3].set_ylim([1e3,1e-15])
    
    for i,ax in enumerate(axs3.flat):
        ax.set_xlabel("He/H number density ratio at atmospheric base", fontsize=18)
    axs3[0].set_ylabel(r"$\rm Mass \;loss\; rate\;[g\;s^{-1}]$", fontsize=18)
    axs3[1].set_ylabel(r"$\rm Fractionation \;factor$", fontsize=18)
    
    masslosses_H = []
    masslosses_He =[]
    
    #
            # For all iontypes and all ions per iontype, find the index in the all_species_data array
            #
            # e.g. iontypes with two types, for H and He would look like
            # iontypes = [     [["H0", 0] , ["p+", 0 ]],
            #                  [["He", 0] , ["He+",0] ],  
            #            ]
    if 0==0:
            for t, iontype in enumerate(iontypes): 
                for s,spc in enumerate(iontype):
                    spc[1] = float('NaN')
            
            for t, iontype in enumerate(iontypes): 
                for s,spc in enumerate(iontype):    #Every iontype is a list of [name, 0], where 0 needs to be replaced by the correct index in the species array
                    for ss,sspc in enumerate(species):  #Every sspc in species is a list [name, mass, charge, namestring]
                        if sspc[0] == spc[0]:
                            #index found
                            spc[1] = ss
            
            for t, iontype in enumerate(iontypes): 
                for s,spc in enumerate(iontype):
                    if math.isnan(spc[1]) == 1:
                        print("ERROR: Can't find species " + spc[0] + " index.")
    
    for i,sim in enumerate(sims):
            timestep = timesteps[i]
            
            timestep_mo_exists = 0
            timestep_po_exists = 0
            
            max_rr  = get_filelength(dirr, "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat" ) #max_rrs[i]
            data_h0 = np.loadtxt(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            if rplanet < 0:
                rp = data_h0[2,0]
            else:
                rp = rplanet
            
            if os.path.isfile(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep-1)+ ".dat"):
                timestep_mo_exists = 1
            if os.path.isfile(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep+1)+ ".dat"):
                timestep_po_exists = 1
            
            defaultr    = 0. * data_h0[:,0]        #An array of 0 elements with the correct radial extent, hence detfault"r"
            defaultdata = 0. * np.asarray(data_h0)
            
            #
            # Setup global array of all species data
            #
            all_species_datas = []
            for s, spc in enumerate(species):
                file = "output_"+ sim +"_"+species[s][0]+"_t" + repr(timestep)+ ".dat"
                #print( "loading speciesdata = " + repr( file ))
                
                tmpdata = []
                if os.path.isfile(dirr + file):
                    tmpdata               = np.loadtxt(dirr + file,   skiprows=2, max_rows=max_rr)
                else:
                    tmpdata = 0. * np.loadtxt(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
                
                #all_species_datas = np.append(all_species_datas, tmpdata)
                if s == 0:
                    all_species_datas = [tmpdata]
                else:
                    #all_species_datas = np.append(all_species_datas, tmpdata)
                    all_species_datas = all_species_datas + [tmpdata]
            
            #
            # Set x-axis variable
            #
            r  = all_species_datas[0][:,0]/rp
            if pressure_on_xaxis == 1:
                r  = total_p/1e6
            dr = np.diff(r)*rp
            
            total_n               = defaultr
            total_p               = defaultr
            total_dens_neutrals   = defaultr
            total_dens_ions       = defaultr
            total_charge          = defaultr
            avg_v_ions            = defaultr
            avg_v_neutrals        = defaultr
            avg_T                 = defaultr
            
            #Helium gets special treatment in this notebook
            total_p_heliums       = defaultr
            mmw                   = defaultr
            scaleH                = defaultr
            grav                  = defaultr
            
            heat_fractional_lineint      = [] #May21
            massloss_predictions         = []
            heat_fractional_volint       = []
            cool_fractional_volint       = []
            
            #print("Done with global array setup")
            #print("all_species_datas = " + repr(all_species_datas))
            #print("all_species_datas[0] = " + repr(all_species_datas[0]))
            #print(" len(allspecies)" + repr(len(all_species_datas)))
            #print("all_species_datas[0][0,0] = " + repr(all_species_datas[0][0,0]))
            
            #
            # Compute total and averaged quantities over all species
            #
            for s, spc in enumerate(all_species_datas):
                total_n      = np.add(total_n, spc[:,4])
                total_charge = np.add(total_charge, species[s][2] * spc[:,4]) 
                total_p      = np.add(total_p, spc[:,10])
                #mmw          = np.add(mmw, species[s][1] * spc[:,4])
                mmw          = np.add(mmw, spc[:,1]/amu)
                avg_T        = np.add(avg_T, spc[:,4] * spc[:,12])
                
                if( species[s][2]== 0 ):
                    avg_v_neutrals      = np.add(avg_v_neutrals, spc[:,2])    
                    total_dens_neutrals = np.add(total_dens_neutrals, spc[:,1])
                else:
                    avg_v_ions     = np.add(avg_v_ions, spc[:,2])
                    total_dens_ions= np.add(total_dens_ions, spc[:,1])
                    
                if np.abs(1. - species[s][1]/4.) < 0.1: #Found helium species
                    total_p_heliums = np.add(total_p_heliums, spc[:,10])
                    
                vol = 4.*3.141592*dr*(r[0:-1]*rp)**2
                dr_fluxes = np.multiply(dr, spc[0:-1,22])
                
                tmp_flux = np.sum( dr_fluxes )
                tmp_heat = np.sum(np.multiply(vol, spc[0:-1,22]) )
                tmp_cool = np.sum(np.multiply(vol, spc[0:-1,21]) )
                
                
                #mdots = [get_energylimit_massloss_onlyplanetradius(F=drflux, mplanet=8., rplanet=r[i]*dr) for i,drflux in enumerate(dr_fluxes)  ]
                
                mdots = np.sum(np.multiply(vol, spc[0:-1,22]) )
                
                for k,rad in enumerate(r):
                    if rad > 1e10: #approximately cs
                        mdots[k] = 0.
                
                tmp_mdot = np.sum(mdots)/(6.678e-8*8.*6e27 / (r[0]*rp) )
                
                heat_fractional_lineint = np.append(heat_fractional_lineint, tmp_flux)
                massloss_predictions    = np.append(massloss_predictions, tmp_mdot)
                
                heat_fractional_volint  = np.append(heat_fractional_volint, tmp_heat)
                cool_fractional_volint  = np.append(cool_fractional_volint, tmp_cool)
            
            with np.printoptions(threshold=np.inf):
                for s,spc in enumerate(species):
                    #print("s / heat lineint fluxes = " + repr([spc[0], heat_fractional_lineint[s] ]))
                    print("s / mdot = " + repr(spc[0]) + f"{Decimal(massloss_predictions[s]):.2E}" )
            
            charge_ratio    = 1. + np.divide(total_charge , total_n)
            avg_v_neutrals /= total_dens_neutrals
            avg_v_ions     /= total_dens_ions
            
            avg_T          /= total_n
            mmw            /= total_n
            
            grav   = -all_species_datas[0][:,19] / (all_species_datas[0][:,0])
            scaleH = kb * avg_T / (mmw * amu * grav)
            
            
            #print("After the index finding routine, iontypes = " + repr(iontypes))
            #
            # Compute total quantities for all iontypes
            #
            timestep_mo_exists = 0
            iontype_data = [] # Will be set of [ion_fraction, mixing_ratio, mass_fluxes] per iontype
            for t, iontype in enumerate(iontypes):
                t_numdens                 = defaultr
                t_massdens                = defaultr
                t_actualions_numfraction  = defaultr
                t_actualions_massfraction = defaultr
                
                t_massflux = defaultr
                t_massflux_po = defaultr
                t_massflux_mo = defaultr
                
                for s,spc in enumerate(iontype):
                    t_massdens    = np.add(t_massdens, all_species_datas[ spc[1] ][:,1]  )
                    t_numdens     = np.add(t_numdens,  all_species_datas[ spc[1] ][:,4]  )
                    t_massflux    = np.add(t_massflux, all_species_datas[ spc[1] ][:,5]  )
                    
                    if( species[ spc[1] ][2] != 0 ): #If charge [2] in species index spc[1] is nonzero
                        t_actualions_massfraction = np.add(t_actualions_massfraction, all_species_datas[ spc[1] ][:,1])
                        t_actualions_numfraction  = np.add(t_actualions_numfraction,  all_species_datas[ spc[1] ][:,4])
                    
                    if(timestep_po_exists):
                        tmp_data        = np.loadtxt(dirr + "output_"+ sim +"_"+spc[0]+"_t" + repr(timestep+1)+ ".dat",   skiprows=2, max_rows=max_rr)
                        t_massflux_po   = np.add(t_massflux_po, tmp_data[:,5])
                    if(timestep_mo_exists):
                        tmp_data        = np.loadtxt(dirr + "output_"+ sim +"_"+spc[0]+"_t" + repr(timestep-1)+ ".dat",   skiprows=2, max_rows=max_rr)
                        t_massflux_mo   = np.add(t_massflux_mo, tmp_data[:,5])
                
                t_ionfrac         = all_species_datas[ iontype[0][1] ][:,4] / t_numdens
                t_mixingratio     = all_species_datas[ iontype[0][1] ][:,4] / total_n
                t_mixingratio_tot = t_numdens / total_n
                t_ionmassratio    = t_actualions_massfraction / t_massdens
                
                
                if t == 0:
                    iontype_data =                 [[t_ionfrac, t_mixingratio, t_massflux, t_massflux_mo, t_massflux_po, t_mixingratio_tot, t_massdens, t_ionmassratio] ]
                else:
                    iontype_data = iontype_data +  [[t_ionfrac, t_mixingratio, t_massflux, t_massflux_mo, t_massflux_po, t_mixingratio_tot, t_massdens, t_ionmassratio] ]
            
            #print("len(iontype data) =" +  repr(len(iontype_data)))
            #print("iontype+data = " + repr(iontype_data))
            #FEB18th
            #
            # Start plotting
            #
            axs0[0,0].set_ylabel(r"$\rm Species \; ratio\;\;  n_s/\sum_{group} n_s$",fontsize=18)
            axs0[1,0].set_ylabel(r"$\rm Mixing\; ratio\;\;  n_s/\sum_{all} n_s$",fontsize=18)
            axs0[2,0].set_ylabel("Mass-loss rate [g/s]",fontsize=18)
            axs0[3,0].set_ylabel("Velocities [cm/s]",fontsize=18)
            
            fracs = []
            
            groupstyles = ['-','--',':','-.']
            for t,iontype in enumerate(iontypes):
                axs0[0,t].set_title(ionnames_total[t],fontsize=18)
                axs0[0,t].loglog(r, iontype_data[t][0], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])      #species fraction among total ion
                axs0[1,t].loglog(r, iontype_data[t][1], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])      #species mixing ratio
                axs0[1,t].loglog(r, iontype_data[t][5], lw=1., c=cols[i], ls=':')            #total mixing ratio
                axs0[2,t].loglog(r, iontype_data[t][2], lw=3., c=cols[i], ls=sstyls[i], label=labels[i])      #mass flux
                
                if(timestep_mo_exists):
                    axs0[2,t].fill_between(r, iontype_data[t][2], iontype_data[t][3], color=cols[i], alpha=0.2, hatch='-') 
                if(timestep_po_exists):
                    axs0[2,t].fill_between(r, iontype_data[t][2], iontype_data[t][4], color=cols[i], alpha=0.2, hatch='-') 
                
                fluxratio_at_bound =  iontype_data[t][2][-1] /  iontype_data[0][2][-1]
                densratio_at_bound =  iontype_data[t][6][2]  /  iontype_data[0][6][0]
                
                if t==2:
                    masslosses_H = np.append(masslosses_H,   iontype_data[t][2][-1])
                if t==3:
                    masslosses_He = np.append(masslosses_He, iontype_data[t][2][-1])
                
                #print("H fractionation factor: " + repr(fluxratio_at_bound/densratio_at_bound) + " fluxratio, densratio = " + repr([fluxratio_at_bound, densratio_at_bound]) )
                fracs = np.append(fracs, fluxratio_at_bound/densratio_at_bound)
                
                for s,spc in enumerate(iontype):
                    #axs0[3,t].loglog(r,     all_species_datas[ spc[1] ][:, 11], c=cols[i], ls=groupstyles[s], lw=2., label= spc if i==0 else ""  )
                    axs0[3,t].loglog(r,     all_species_datas[ spc[1] ][:, 11], c=cols[i], ls=groupstyles[s], lw=2., label= species[spc[1]][3] if i==0 else ""  )
                
            print("All fractionation factors for groups: " + repr(fracs))
            #print("All H ")
            
            axs1[0,0].loglog(r, avg_T, lw=2., c=cols[i], ls=sstyls[i], label=labels[i])      #C/O ratio in flow
            axs1[0,1].loglog(r, avg_v_neutrals, lw=2., c=cols[i], ls=sstyls[i], label=labels[i])
            axs1[0,1].loglog(r, avg_v_ions, lw=1., c=cols[i], ls=sstyls[i], label='')
            axs1[0,2].loglog(r, total_p/1e6,         lw=3., c=cols[i], ls=sstyls[i], label=labels[i])     
            axs1[0,2].loglog(r, total_p_heliums/1e6,         lw=3., c=cols[i], ls=':')     
            
            densratio    =  iontype_data[plot_group_massloss_ratio[0]][6][0]/iontype_data[plot_group_massloss_ratio[1]][6][0] #For group 0 find iontype data, extract column 6 (mass density) and take cell 2, repeat for group 1, divide
            numdensratio =  iontype_data[plot_group_massloss_ratio[0]][6][0]/iontype_data[plot_group_massloss_ratio[1]][6][0] * 0.25
            mdotratio    =  iontype_data[plot_group_massloss_ratio[0]][2]   /iontype_data[plot_group_massloss_ratio[1]][2]
            mdotratio2   =  iontype_data[plot_group_massloss_ratio[0]][2][-2]/iontype_data[plot_group_massloss_ratio[1]][2][-2]
            
            axs1[0,3].loglog(r, mdotratio / densratio,             lw=3, c=cols[i], ls=sstyls[i], label=labels[i])  
            axs1[1,0].loglog(r, scaleH / rp, lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            axs1[1,1].semilogx(r, mmw, lw=3, c=cols[i], ls=sstyls[i], label=labels[i])    
            axs1[1,2].loglog(r, grav, lw=3, c=cols[i], ls=sstyls[i], label=labels[i]) 
            
            #axs1[1,3].loglog(r, charge_ratio, lw=1.5, c=cols[i], ls=sstyls[i], label=labels[i])
            axs1[1,3].loglog(avg_T, total_p/1e6, lw=2.5, c=cols[i], ls=sstyls[i], label=labels[i])
            
            newlist = [elm if elm > 1e-1 else 1e-1 for elm in avg_v_neutrals]
            #ax77.loglog(newlist, total_p/1e6, lw=2.5, c=cols[i], ls='--', label=labels[i])
            ax77.loglog(newlist, total_p/1e6, lw=2.5, c=cols[i], ls='', label=labels[i], marker='o', markersize='7')
            
            #for s, spc in enumerate(all_species_datas):
            for sn,snap in enumerate(snapshot):
                if i == snap:
                    #for s, spc in enumerate(species):
                    #    file = "output_"+ sim +"_"+species[s][0]+"_t" + repr(timestep)+ ".dat"

                    neutrcnt = 0
                    chargecnt= 0
                    for s, spc in enumerate(species):
                        if spc[2] == 0:
                            axs2[sn,0].loglog(r, all_species_datas[s][:,4] / total_n, label = spc[3], lw=3., ls= '-' if neutrcnt%2==0 else '--' )
                            neutrcnt += 1
                        else:
                            axs2[sn,1].loglog(r, all_species_datas[s][:,4] / total_n, label = spc[3], lw=3., ls= '--' if chargecnt > 9 else '-' )
                            chargecnt += 1
            
            
            for a, ax in enumerate(axs2[:,0]):
                axs2[a,0].set_title(labels[snapshot[a]], fontsize=13)
                axs2[a,1].set_title(labels[snapshot[a]], fontsize=13)
                
                axs2[a,0].legend(fontsize=13, ncol=2, bbox_to_anchor=(1.0, 1.00), loc="upper left") #APRIL3
                axs2[a,1].legend(fontsize=13, ncol=2, bbox_to_anchor=(1.0, 1.00), loc="upper left")
                
                #add pressure levels
                #pressures = [1., 1e-3,1e-6,1e-9,1e-12]
                #print(" isobars mask " + repr(plot_isobars_mask[a]) + " len " + repr(len(plot_isobars_mask[a])))
                #print(" isobars labels " + repr(plot_isobars_labels) + " len " + repr(len(plot_isobars_labels)))
                if len(plot_isobars_mask[a]) > 0:
                    if i == snapshot[a]:
                        ff = 0.999
                        if snapshot[a] >= len(sims)-2:
                            pressures = plot_isobars_mask[a] #[1e-3, 1e-6, 1e-9]
                            ff=0.999
                        else:
                            pressures = plot_isobars_mask[a] #[1e-3, 1e-6, 1e-9]
                        
                        plabels   = ["0.1 mbar", r"$\rm \mu bar$", "nbar" , "pbar"]
                        if len(plot_isobars_labels) > 0:
                            plabels   = plot_isobars_labels[a]
                            
                        for p,pressure in enumerate(pressures):

                            iso_r = r[np.where(total_p/1e6 > pressure)]
                            #print(" looking for pressure = " + repr(pressure) + " with all_ps = " + repr(total_p/1e6))
                            iso_r = iso_r[-1]
                            axs2[a,0].axvline(x=iso_r, c='lightgrey', ls='--')
                            axs2[a,1].axvline(x=iso_r, c='lightgrey', ls='--')
                            axs2[a,0].text(s=plabels[p], x=iso_r*ff, y=1e-3, c='darkgrey', fontsize=14, rotation=90)
                            axs2[a,1].text(s=plabels[p], x=iso_r*ff, y=1e-3, c='darkgrey', fontsize=14, rotation=90)
            
            for t,iontype in enumerate(iontypes):
                mdot_total  = iontype_data[t][2][-2]
                ionfraction = iontype_data[t][7][-2]
                name        = ionnames_total[t]
                
                axs3[0].loglog(numdensratio, mdot_total,             c=cols[t], ls="", marker="o", markersize="16", label=name+" all " if i==0 else "") #Total mass loss for this group at edge, replace -2 with sonic point
                axs3[0].loglog(numdensratio, mdot_total*ionfraction, c=cols[t], ls="", marker="+", markersize="16", label=name+" ion " if i==0 else "") #Total mass loss for this group at edge, replace -2 with sonic point
            axs3[1].loglog(numdensratio,     mdotratio2 / densratio, lw=3, c=cols[i], ls="", marker="^", markersize="16", label=labels[i])    
            #axs3[1].loglog(numdensratio,     mdotratio2 , lw=3, c=cols[i], ls="", marker="+", markersize="16")    
            #axs3[1].loglog(numdensratio,     densratio, lw=3, c=cols[i], ls="", marker="o", markersize="16")    
            
            for a, ax in enumerate(axs0.flat):
                ax.legend(fontsize=11)
            for a, ax in enumerate(axs0[3,:]):
                ax.legend(fontsize=14)
            for a, ax in enumerate(axs1.flat):
                ax.legend(fontsize=11)
            for a, ax in enumerate(axs2.flat):
                ax.legend(fontsize=13, ncol=2)
            for a, ax in enumerate(axs2[:,0]):
                #ax.set_title(labels[a], fontsize=13)
                #ax.set_title(labels[snapshot[a]], fontsize=13)
                axs2[a,1].set_title(labels[snapshot[a]], fontsize=13)
    
    print(" All mdot H = " + repr(masslosses_H))
    print(" All mdot He = " + repr(masslosses_He))
    
    print("Writing he_frac_masslossrates.txt...")
    filehandle = open("he_frac_masslossrates.txt", "a")
    filehandle.writelines("#"+outlabel+"\n")
    #filebuffer = ['%i %1.8e %1.8e \n'%(elm, masslosses_H[i], masslosses_He[i]) for i,elm in enumerate(pradii) ]
    #filehandle.writelines(filebuffer)
    filehandle.writelines(" ")
    filehandle.close()
    
    xlabelstring = r"$\rm Radius \; [R_p]$"
    if rplanet > 0: #Assuming the radius given is in earth radii
        xlabelstring = r"$\rm Radius \; [R_{\oplus}]$"
        axs1[1,0].set_ylabel(r"$\rm Scale\;Height \; [R_{\oplus}]$", fontsize=16)
    if pressure_on_xaxis == 1:
        xlabelstring = r"$\rm Pressure \; [bar]$"
    
    for a, ax in enumerate(axs0[3,:]):
                ax.set_xlabel(xlabelstring,fontsize=16)
    for a, ax in enumerate(axs1[0,:]):
                ax.set_xlabel(xlabelstring,fontsize=16)
    for a, ax in enumerate(axs1[1,:]):
                ax.set_xlabel(xlabelstring,fontsize=16)
    for a, ax in enumerate(axs2[1,:]):
                ax.set_xlabel(xlabelstring,fontsize=16)
    for a, ax in enumerate(axs3.flat):
                ax.legend(fontsize=13, ncol=2)
            
    axs1[1,3].set_xlabel("Temperature [K](-), Velocity [cm/s](o)",fontsize=16)
    axs1[1,3].set_xlim([450,20000.])
    axs1[1,3].set_ylim([1e0,1e-12])
    ax77.set_xlim([1e0,1e5])

    
    #axs2.legend(fontsize=11)        
    #axs[1,2].axhline(y=1,c='k')
    #spl5.axhline(y=1,c='k')
    #FEB4
    
    fig0.tight_layout()
    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    
    return fig0, fig1, fig2, fig3


#
# just prints the mass loss data for species groups g0, g1, g2 .. gfinal, on a multidimensional grid for varying parameters p0, p1, p2, .. pfinal etc. in the format
# 
# p0, p1, p2 ,... , pfinal, mdot_g0, mdot_g1, mdot_g2, ... mdot_gfinal
#
def print_mdots_for_groups(directory, sims, times, labels, paramsets, styles, title,  
                                           species, iontypes, ionnames_total ):
    
    num_ions = len(iontypes)
    if(num_ions < 1):
        print("ERROR: num_ions has to be >= 1!   num_ions = " + repr(num_ions))
    #else:
    #    print(" num_ions = " + repr(num_ions))
    
    masslosses_H = []
    masslosses_He =[]
    
    #
            # For all iontypes and all ions per iontype, find the index in the all_species_data array
            #
            # e.g. iontypes with two types, for H and He would look like
            # iontypes = [     [["H0", 0] , ["p+", 0 ]],
            #                  [["He", 0] , ["He+",0] ],  
            #            ]
    if 0==0:
            for t, iontype in enumerate(iontypes): 
                for s,spc in enumerate(iontype):
                    spc[1] = float('NaN')
            
            for t, iontype in enumerate(iontypes): 
                for s,spc in enumerate(iontype):    #Every iontype is a list of [name, 0], where 0 needs to be replaced by the correct index in the species array
                    for ss,sspc in enumerate(species):  #Every sspc in species is a list [name, mass, charge, namestring]
                        if sspc[0] == spc[0]:
                            #index found
                            spc[1] = ss
            
            for t, iontype in enumerate(iontypes): 
                for s,spc in enumerate(iontype):
                    if math.isnan(spc[1]) == 1:
                        print("ERROR: Can't find species " + spc[0] + " index.")
    
    for i,sim in enumerate(sims):
            timestep = times[i]
            
            timestep_mo_exists = 0
            timestep_po_exists = 0
            
            max_rr  = get_filelength(dirr, "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat" ) #max_rrs[i]
            data_h0 = np.loadtxt(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
            
            rp = data_h0[2,0]
            
            defaultr    = 0. * data_h0[:,0]        #An array of 0 elements with the correct radial extent, hence detfault"r"
            defaultdata = 0. * np.asarray(data_h0)
            
            #
            # Setup global array of all species data
            #
            all_species_datas = []
            for s, spc in enumerate(species):
                file = "output_"+ sim +"_"+species[s][0]+"_t" + repr(timestep)+ ".dat"
                #print( "loading speciesdata = " + repr( file ))
                
                tmpdata = []
                if os.path.isfile(dirr + file):
                    tmpdata               = np.loadtxt(dirr + file,   skiprows=2, max_rows=max_rr)
                else:
                    tmpdata = 0. * np.loadtxt(dirr + "output_"+ sim +"_"+species[0][0]+"_t" + repr(timestep)+ ".dat",   skiprows=2, max_rows=max_rr)
                
                #all_species_datas = np.append(all_species_datas, tmpdata)
                if s == 0:
                    all_species_datas = [tmpdata]
                else:
                    #all_species_datas = np.append(all_species_datas, tmpdata)
                    all_species_datas = all_species_datas + [tmpdata]
            
            #
            # Set x-axis variable
            #
            r  = all_species_datas[0][:,0]/rp
            dr = np.diff(r)*rp
            
            total_n               = defaultr
            total_p               = defaultr
            total_dens_neutrals   = defaultr
            total_dens_ions       = defaultr
            total_charge          = defaultr
            avg_v_ions            = defaultr
            avg_v_neutrals        = defaultr
            avg_T                 = defaultr
            
            #Helium gets special treatment in this notebook
            total_p_heliums       = defaultr
            mmw                   = defaultr
            scaleH                = defaultr
            grav                  = defaultr
            
            heat_fractional_lineint      = [] #May21
            massloss_predictions         = []
            heat_fractional_volint       = []
            cool_fractional_volint       = []
            
            #print("Done with global array setup")
            #print("all_species_datas = " + repr(all_species_datas))
            #print("all_species_datas[0] = " + repr(all_species_datas[0]))
            #print(" len(allspecies)" + repr(len(all_species_datas)))
            #print("all_species_datas[0][0,0] = " + repr(all_species_datas[0][0,0]))
            
            #
            # Compute total and averaged quantities over all species
            #
            for s, spc in enumerate(all_species_datas):
                total_n      = np.add(total_n, spc[:,4])
                total_charge = np.add(total_charge, species[s][2] * spc[:,4]) 
                total_p      = np.add(total_p, spc[:,10])
                #mmw          = np.add(mmw, species[s][1] * spc[:,4])
                mmw          = np.add(mmw, spc[:,1]/amu)
                avg_T        = np.add(avg_T, spc[:,4] * spc[:,12])
                
                if( species[s][2]== 0 ):
                    avg_v_neutrals      = np.add(avg_v_neutrals, spc[:,2])    
                    total_dens_neutrals = np.add(total_dens_neutrals, spc[:,1])
                else:
                    avg_v_ions     = np.add(avg_v_ions, spc[:,2])
                    total_dens_ions= np.add(total_dens_ions, spc[:,1])
                    
                if np.abs(1. - species[s][1]/4.) < 0.1: #Found helium species
                    total_p_heliums = np.add(total_p_heliums, spc[:,10])
                    
                vol = 4.*3.141592*dr*(r[0:-1]*rp)**2
                dr_fluxes = np.multiply(dr, spc[0:-1,22])
                
                tmp_flux = np.sum( dr_fluxes )
                tmp_heat = np.sum(np.multiply(vol, spc[0:-1,22]) )
                tmp_cool = np.sum(np.multiply(vol, spc[0:-1,21]) )
                
                mdots = np.sum(np.multiply(vol, spc[0:-1,22]) )
                
                for k,rad in enumerate(r):
                    if rad > 1e10: #approximately cs
                        mdots[k] = 0.
                
                tmp_mdot = np.sum(mdots)/(6.678e-8*8.*6e27 / (r[0]*rp) )
                
                heat_fractional_lineint = np.append(heat_fractional_lineint, tmp_flux)
                massloss_predictions    = np.append(massloss_predictions, tmp_mdot)
                
                heat_fractional_volint  = np.append(heat_fractional_volint, tmp_heat)
                cool_fractional_volint  = np.append(cool_fractional_volint, tmp_cool)
            
            #with np.printoptions(threshold=np.inf):
            #    for s,spc in enumerate(species):
            #        #print("s / heat lineint fluxes = " + repr([spc[0], heat_fractional_lineint[s] ]))
            #        #print("s / mdot = " + repr(spc[0]) +" "+ f"{Decimal(massloss_predictions[s]):.2E}" )
            
            charge_ratio    = 1. + np.divide(total_charge , total_n)
            avg_v_neutrals /= total_dens_neutrals
            avg_v_ions     /= total_dens_ions
            
            avg_T          /= total_n
            mmw            /= total_n
            
            grav   = -all_species_datas[0][:,19] / (all_species_datas[0][:,0])
            scaleH = kb * avg_T / (mmw * amu * grav)
            
            #print("After the index finding routine, iontypes = " + repr(iontypes))
            #
            # Compute total quantities for all iontypes
            #
            timestep_mo_exists = 0
            iontype_data = [] # Will be set of [ion_fraction, mixing_ratio, mass_fluxes] per iontype
            for t, iontype in enumerate(iontypes):
                t_numdens                 = defaultr
                t_massdens                = defaultr
                t_actualions_numfraction  = defaultr
                t_actualions_massfraction = defaultr
                
                t_massflux = defaultr
                t_massflux_po = defaultr
                t_massflux_mo = defaultr
                
                for s,spc in enumerate(iontype):
                    t_massdens    = np.add(t_massdens, all_species_datas[ spc[1] ][:,1]  )
                    t_numdens     = np.add(t_numdens,  all_species_datas[ spc[1] ][:,4]  )
                    t_massflux    = np.add(t_massflux, all_species_datas[ spc[1] ][:,5]  )
                    
                    if( species[ spc[1] ][2] != 0 ): #If charge [2] in species index spc[1] is nonzero
                        t_actualions_massfraction = np.add(t_actualions_massfraction, all_species_datas[ spc[1] ][:,1])
                        t_actualions_numfraction  = np.add(t_actualions_numfraction,  all_species_datas[ spc[1] ][:,4])
                    
                    if(timestep_po_exists):
                        tmp_data        = np.loadtxt(dirr + "output_"+ sim +"_"+spc[0]+"_t" + repr(timestep+1)+ ".dat",   skiprows=2, max_rows=max_rr)
                        t_massflux_po   = np.add(t_massflux_po, tmp_data[:,5])
                    if(timestep_mo_exists):
                        tmp_data        = np.loadtxt(dirr + "output_"+ sim +"_"+spc[0]+"_t" + repr(timestep-1)+ ".dat",   skiprows=2, max_rows=max_rr)
                        t_massflux_mo   = np.add(t_massflux_mo, tmp_data[:,5])
                
                t_ionfrac         = all_species_datas[ iontype[0][1] ][:,4] / t_numdens
                t_mixingratio     = all_species_datas[ iontype[0][1] ][:,4] / total_n
                t_mixingratio_tot = t_numdens / total_n
                t_ionmassratio    = t_actualions_massfraction / t_massdens
                
                if t == 0:
                    iontype_data =                 [[t_ionfrac, t_mixingratio, t_massflux, t_massflux_mo, t_massflux_po, t_mixingratio_tot, t_massdens, t_ionmassratio] ]
                else:
                    iontype_data = iontype_data +  [[t_ionfrac, t_mixingratio, t_massflux, t_massflux_mo, t_massflux_po, t_mixingratio_tot, t_massdens, t_ionmassratio] ]
            
            fracs = []
            
            groupstyles = ['-','--',':','-.']
            #grpstr = [repr(iontype_data[g][2][-1])  for g,group in enumerate(iontypes) ]
            #groupstr = [ iontype_data[g][2][-1]  for g,group in enumerate(iontypes) ] #JUNE25
            groupstr = [f'{Decimal(iontype_data[g][2][-1]):.2E}'  for g,group in enumerate(iontypes) ]
            
            #print(repr(paramsets[i]) + " " + repr(groupstr))
            
            #print(' '.join([str(param) for p,param in enumerate(paramsets[i])   ]) + " " + ' '.join(groupstr))
            print(' '.join([str(param) if p<3 else f'{Decimal(param):.2E}' for p,param in enumerate(paramsets[i])   ]) + " " + ' '.join(groupstr))
            
            
#
#
# Start printing
#
#
#


rp =1.37e8
title=''
cols = [ 'red', 'blue', 'green', 'orange', 'magenta', 'cyan', 'gold', 'olive', 'lightgrey']
sstyls = ['-','-',     '-',      '-',    '-',       '-',         '-',  '-',     '-']

dirr = "/home/mschulik/aiolos/runs_2023_co/"
labels  = []

#Species array is an array of [name, mass, charge], name must be the same as in the simulation
species = [["S0", 1., 0, r"$\rm H^0$"],     ["S1", 1, +1, r"$\rm H^+$"],       ["e-", 5e-4, -1, r"$\rm e^-$"],
           ["He11S", 4., 0, r"$\rm He^0$"], ["He23S", 4, 0, r"$\rm He2^3S$"],  ["Hep", 4., +1,r"$\rm He^+$"],  ["Hepp", 4., +2,r"$\rm He^{++}$"] ,
            ["H2", 2,  0, r"$\rm H_2$"],  ["H2p", 2, +1, r"$\rm H_2^{+}$"],
          ]
#Iontypes groups ions and neutrals that belong together, in order to plot neutral or ion fractions. The 0 index is a dummy and needs to be there
groups = [
              [["H2", 0], ["H2p", 0]  , ["S0", 0] , ["S1", 0 ]], 
              [["He23S", 0], ["He11S", 0] , ["Hep", 0 ], ["Hepp", 0]],   
           ] 

groupnames_total = [ r"$\rm H_2/all \;H$", r"$\rm He^{+}/all \;He$", "He triplet", r"$\rm He^{0}/(H_2 + H^{0})$"]
outlabel = ""
#
# Finish setup for all
#
s1 = "8me-quicktest_G2-10days-r"
s2 = "-noshad-H2-long-He"
s3 = "-kzz"

gmplanet = 6.678e-8*8.*5.98e27

#
# Nominal FUV flux
#
timesteps = [[5,  4,-1, 4,  4, 4, -1, -1],
             [-1,-1,-1, -1,-1,-1, -1,-1],
             [-1,-1,-1, -1,-1,-1, -1,-1]]

for r,r00 in enumerate([20,30,40]):
    r0  = r00
    fuv = 120
    heliums    = [40, 60, 80, 90, 95, 99, 999, 9999]
    heliums100 = [0.40, 0.60, 0.80, 0.90, 0.95, 0.99, 0.999, 0.9999]
    hefracs = [1./(1./he-1.) for he in heliums100]
    T     = 1000.
    
    pbar_to_pcgs = 1e6
    pmix  = 1e-6 * pbar_to_pcgs #pmix = 1 microbar

    ncrit = pmix/(kb*T)
    nH2   = [(1-he)*ncrit for he in heliums100]
    nHe   = [he*ncrit for he in heliums100]
    b     = 5e17*T**0.75
    kzz   = b/ncrit

    mdot_dl = [(4.-2.)*amu*gmplanet*b*(1.-he)/(kb*T)*(2.*amu)   for he in heliums100  ]

    paramsets = [[r0, heliums[h], hefracs[h], fuv, kzz, mdot_dl[h]] for h,he in enumerate(hefracs)]
    sims = [s1+repr(pset[0])+s2+repr(pset[1])+s3 for pset in paramsets]
    
    print_mdots_for_groups(dirr, sims, timesteps[r], labels, paramsets, sstyls, title, species, groups, groupnames_total)

#
# x10 FUV flux
#
s3 = "-kzz-fuv10"
timesteps = [[1, 1,5,  1,1,1, -1, -1],
             [0, 1,1,  2,2,2, -1,-1],
             [2, -1,-1,-1,-1,-1 ,-1,-1],
             [-1,-1,-1, -1,-1,-1, -1,-1], ]
for r,r00 in enumerate([20,30,40,50]):
    r0  = r00
    fuv = 1200
    heliums    = [40, 60, 80, 90, 95, 99, 999, 9999]
    heliums100 = [0.40, 0.60, 0.80, 0.90, 0.95, 0.99, 0.999, 0.9999]
    hefracs = [1./(1./he-1.) for he in heliums100]
    T     = 1000.

    pbar_to_pcgs = 1e6 
    pmix  = 1e-6 * pbar_to_pcgs #pmix = 1 microbar
    ncrit = pmix/(kb*T)
    nH2   = [(1-he)*ncrit for he in heliums100]
    nHe   = [he*ncrit for he in heliums100]
    b     = 5e17*T**0.75
    kzz   = b/ncrit
    
    mdot_dl = [(4.-2.)*amu*gmplanet*b*(1.-he)/(kb*T)*(2.*amu)   for he in heliums100  ]
    
    paramsets = [[r0, heliums[h], hefracs[h], fuv, kzz, mdot_dl[h]] for h,he in enumerate(hefracs)]
    sims = [s1+repr(pset[0])+s2+repr(pset[1])+s3 for pset in paramsets]
    
    print_mdots_for_groups(dirr, sims, timesteps[r], labels, paramsets, sstyls, title, species, groups, groupnames_total)
