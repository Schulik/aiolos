import numpy as np

#plot for [i=-1:-1:1] "runs_guillot2010/output_hd209458b_nominal_1s1b_H0_t".i.".dat" u 1:13 w l lw 2 t "nominal, o ".i, 
#for [i=-1:-1:2] "runs_guillot2010/output_hd209458b_twoconstant_1s1b_H0_t".i.".dat" u 1:13 w l lw 2 t "y=10., o ".i, 
#for [i=-1:-1:1] "runs_guillot2010/output_hd209458b_twoconstant2_1s1b_H0_t".i.".dat" u 1:13 w l lw 2 t "y=0.25, o ".i, 
#for [i=-1:-1:1] "runs_guillot2010/output_hd209458b_twoconstantRout_1s1b_H0_t".i.".dat" u 1:13 w l lw 2 t "y=10, large o ".i, 
#for [i=-1:-1:1] "runs_guillot2010/output_hd209458b_twoconstant2Rout_1s1b_H0_t".i.".dat" u 1:13 w l lw 2 t "y=0.25, large o ".i, 
#"runs_guillot2010/output_hd209458b_twoconstant_1s2b_H0_t-1.dat" u 1:13 w l lw 2 t "y=10, 2band", 
#"runs_guillot2010/output_hd209458b_twoconstant_1s2b_l06_H0_t2.dat" u 1:13 w l lw 2 t "y=10., 2band l06", 
#"runs_guillot2010/output_hd209458b_twoconstant_1s2b_l20_H0_t-1.dat" u 1:13 w l lw 2 t "y=10., 2band, l20"

# Analytic comparison w/ guillot
# plot "runs_guillot2010/diagnostic_hd209458b_twoconstantRout_1s1b_t-1.dat" u 9:14 w l lw 2 t "y=10, large o -1", "runs_guillot2010/diagnostic_hd209458b_twoconstant2Rout_1s1b_t-1.dat" u 9:14 w l lw 2 t "y=0.25, large o -1", "runs_guillot2010/analytic_g2010_gamma10.dat" u 1:2 w l lw 2 t "analytic, y=10", "runs_guillot2010/analytic_g2010_gamma0.25.dat" u 1:2 w l lw 2 t "analytic, y=0.25"
#
Tirr = 1431.
Tint = 100.

kv   = 1.e-1
mustar = np.cos(0.)
kvstar = kv/mustar
k1   = 1.e-2
k2   = 1.e-2
beta = 0. #...1
R      = k1/k2 

kr     = k1*k2 / (beta*k2 + (1.-beta)*k1)
kp     = beta*k1 + (1.-beta)*k2

gammap = kp/kr
gammav = kv/kr
gammavstar = kvstar/mustar
gamma1 = k1/kr
gamma2 = k2/kr
gamma  = kv/k1
taulim = 1./(gamma1*gamma2)*np.sqrt(gammap/3.)

#a0 = 1./gamma1 + 1./gamma2
#a1 = (gamma1+gamma2)*taulim - (...)
#a1 = -1./(3.*taulim**2.) ()

#Grey limit
A = 2./3.
B = 0
C = 2./3. - 2./gammavstar**2. + 2./gammavstar   + 2.*np.log(1+gammavstar) * (1./gammavstar**3. - 1./(3*gammavstar))
D = 0
E = gammavstar/3. - 1./gammavstar

print("Taulim, tau=")
print(taulim)
print(gammavstar)
print("A,b,c,d,e")
print(A)
print(C)
print(E)

taulist = 10.**np.linspace(-3,4,100)


#temperatures = (3./4.*Tint**4.*(tau + A + B*np.exp(-tau/taulim)) + 3./4.*Tirr**4.*mustar*(C + D*np.exp(-tau/taulim) + E * np.exp(-gammavstar * tau) ))**0.25
#temperatures = [ (3./4.*Tint**4.*(tau + A + B*np.exp(-tau/taulim)) + 3./4.*Tirr**4.*mustar*(C + D*np.exp(-tau/taulim) + E * np.exp(-gammavstar * tau) ))**0.25 for tau in taulist ]
temperatures = [  ( 3./4.*Tint**4.*(tau + A + B*np.exp(-tau/taulim)) + 3./4.*Tirr**4.*mustar*(1.*C +  1.*D*np.exp(-tau/taulim) + 1.* E * np.exp(-gammavstar * tau)  ))**0.25 for tau in taulist ]

sq3 = np.sqrt(3.)

temperatures2 =  [( 3./4.*Tint**4.*(2./3. + tau) + 3./4.*Tirr**4.* 1./4. * (2./3. + 1./(gamma*sq3) + (gamma/sq3 - 1./(gamma*sq3))*np.exp(-gamma*tau*sq3 )) )**0.25 for tau in taulist ]

datas = zip(taulist,temperatures2)


#print(list(datas))
for i,j in list(datas):
    print(i,j)
    
print(gamma)
