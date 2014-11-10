from pylab import *
from scipy import constants as const
import numpy as np
con_h = const.h/const.e
con_k = const.k/const.e

np.seterr(all = "ignore")

parsec = 3.08567758e16

m1 = 5.0 			# central mass
alfa = 0.5			# accretion parameter
mass_flow = 5.0		# mass flow from disk
Rc = 1.0*10.0**-4		# diameter of central object /10**10cm

Ri = 1.001*10.0**-4		# inner radius
Rf = 10.0**-3			# outer radius
Rstep = 10.0**-6		# step for computing

Tc = []
R_x = []

for R in np.arange(Ri,Rf,Rstep):		# creating the temperature profile
	R_x.append(R)

	f = (1.0-(Rc/R)**(0.5))**(0.25)
	Tc.append(1.4*10.0**4.0*alfa**(-1.0/5)*mass_flow**(3.0/10)*m1**(1.0/4)*R**(-3.0/4)*f**(6.0/5))

figure(figsize=(8,5))
title('T$_{c}$')
xlabel('R$_{10}$ (R/($10^{10}$ cm))')
ylabel('T$_{c}$ (K)')
plot(R_x,Tc)
savefig("Tc.png")

def planck_law(E,T_in,T_out,step = 1):			#returns rectangular integral of the planck_law created for constant energy and changing temperature 
	
	
	c1 = 2.0*const.pi/(con_h**3.0*const.c**2.0)*const.e

	total = 0.0
	planck = []
	x = []
	
	i = float(T_out)
	while i<=float(T_in):
		B = (i/T_in)**(-11.0/3.0)/T_in*c1*(E)**3.0*1.0/(np.exp(E/(con_k*i))-1.0)
		total += B*step

		i+=step
	
	
	
	return total

def mbb(r_in,E_in,E_out,E_step,T_in,T_out,T_step,D,theta):				# construct the mbb spectrum, takes r_in and T_in from the temperature profile
	
	c1 = 8.0*const.pi*r_in**2.0*np.cos(np.radians(theta))/(3.0*D**2.0)

	f_E = []
	f = []
	
	i = float(E_in)
	while i<=float(E_out):
		f_temp = c1*planck_law(i,T_in,T_out,T_step)
		print i
		f_E.append(i/1000.0)
		f.append(f_temp)
		i += E_step
	return f_E, f

T_in_arg = np.argmax(Tc)
T_in = Tc[T_in_arg]
r_in = R_x[T_in_arg]*10.0**10.0/100.0

f_E,f = mbb(r_in,0.1,20001.0,10.0,T_in,1.0,1000.0,8.0*parsec,10.0)


figure(figsize=(8,5))
title('Planck')
xlabel('E (keV)')
ylabel('B$_{E}$(T) (W.m$^{-2}$.eV$^{-1}$)')
grid()
xlim([-1,20])
xticks(range(0,20,1))
plot(f_E,f)
savefig("planck.png")

