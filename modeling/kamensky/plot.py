from pylab import *
import numpy as np

m1 = 5.0 			# central mass
alfa = 1.0			# accretion parameter
mass_flow = 0.9			# mass flow from disk
Rc = 2.0*10.0**-4		# diameter of central object /10**10cm

Ri = 2.5*10.0**-4		# inner radius
Rf = 10.0**-2			# outer radius
Rstep = 0.00001			# step for computing

R_x = []

SIGMA = []			# surface density integrated from RHO in z direction
H = []				# scaleheight
RHO = []			# density
Tc = []				# temperature
TAU = []			# opacity
NI = []				# viscosity
Vr = []				# radial speed

for R in np.arange(Ri,Rf,Rstep):
	R_x.append(R)
	f = (1.0-(Rc/R)**(0.5))**(0.25)
	SIGMA.append(5.2*alfa**(-4/5)*mass_flow**(7/10)*m1**(1/4)*R**(-3/4)*f**(14/5))
	H.append(1.7*10**8*alfa**(-1/10)*mass_flow**(3/20)*m1**(-3/8)*R**(9/8)*f**(3/5))
	RHO.append(3.1*10**(-8)*alfa**(-7/10)*mass_flow**(11/20)*m1**(5/8)*R**(-15/8)*f**(11/5))
	Tc.append(1.4*10**4*alfa**(-1/5)*mass_flow**(3/10)*m1**(1/4)*R**(-3/4)*f**(6/5))
	TAU.append(190.0*alfa**(-4/5)*mass_flow**(1/5)*f**(4/5))
	NI.append(1.8*10**14*alfa**(4/5)*mass_flow**(3/10)*m1**(-1/4)*R**(3/4)*f**(6/5))
	Vr.append(2.7*10**4*alfa**(4/5)*mass_flow**(3/10)*m1**(-1/4)*R**(-1/4)*f**(-14/5))
	
matplotlib.rcParams.update({'font.size': 12, 'font.family': 'serif'})

figure(figsize=(8,5))
title('$\Sigma$')
xlabel('R$_{10}$ (R/($10^{10}$ cm))')
ylabel('$\Sigma$ (g.cm$^{-2}$)')
plot(R_x,SIGMA)
savefig("SIGMA.png")


figure(figsize=(8,5))
title('H')
xlabel('R$_{10}$ (R/($10^{10}$ cm))')
ylabel('H (cm)')
plot(R_x,H)
savefig("H.png")


figure(figsize=(8,5))
title('$\\rho$')
xlabel('R$_{10}$ (R/($10^{10}$ cm))')
ylabel('$\\rho$ (g.cm$^{-3}$)')
plot(R_x,RHO)
savefig("RHO.png")


figure(figsize=(8,5))
title('T$_{c}$')
xlabel('R$_{10}$ (R/($10^{10}$ cm))')
ylabel('T$_{c}$ (K)')
plot(R_x,Tc)
savefig("Tc.png")


figure(figsize=(8,5))
title('$\\tau$')
xlabel('R$_{10}$ (R/($10^{10}$ cm))')
ylabel('$\\tau$')
plot(R_x,TAU)
savefig("TAU.png")


figure(figsize=(8,5))
title('$\\nu$')
xlabel('R$_{10}$ (R/($10^{10}$ cm))')
ylabel('$\\nu$ (cm$^{2}$s$^{-1}$)')
plot(R_x,NI)
savefig("NU.png")


figure(figsize=(8,5))
title('V$_{r}$')
xlabel('R$_{10}$')
ylabel('v$_{r}$ (cm.s$^{-1}$)')
plot(R_x,Vr)
savefig("vr.png")