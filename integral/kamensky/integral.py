from pylab import *
from astropy.constants import sigma_sb as SB
from scipy import constants as const
import numpy as np
SB = float(SB.value)

def planck_law(T,lam_i,lam_f,step = 0.1):			# units are m^-2.nm^-1.sr^-1
	
	total = SB*T**4.0					# total output from Stephan-Boltzman law
	planck = []
	x = []
	
	for i in arange(lam_i,lam_f,step):
		planck.append(((2.0*const.h*const.c**2.0)/(i*1e-9)**5.0)*(1.0/(exp(const.h*const.c/(i*1e-9*T*const.k))-1))*1e-9)	# *1e-9 for nm^-1
		x.append(i)
	
	return x,planck, total

def integral1(x, y):						# simple stripes method
	integral = 0.0
	h = x[1]-x[0]						# we presume constant step on x axis
	for i in range(1,len(x)):
		integral += h*y[i]
	return integral

def integral2(x, y):						# trapezium method
	integral = 0.0
	h = x[1]-x[0]
	for i in range(1,len(x)):
		integral += (y[i]+y[i-1])*h/2.0
	return integral

def integral3(x, y):						# simpson method
	integral = 0.0
	h = x[1]-x[0]
	
	if len(x)%2 == 0:
		del y[-1]				# simpson method requires odd number of elements
		del x[-1]				# simple sacrifice of the last element is the number is even
								# did not look for better solution, perhaps testing if last element is relevant
	
	integral += h/3.0*(y[0]+y[-1])				# first and last element with coeficient 1
	
	for i in range(1,len(x),2):				# odd elements with coeficient 4
		integral += h/3.0*4.0*y[i]
		

	for i in range(2,len(x)-1,2):				# even elements with coeficient 2
		integral += h/3.0*2.0*y[i]
		
	
	return integral

planck_x, planck, total = planck_law(5000,5,10000,0.01)

int1 = integral1(planck_x,planck)
int2 = integral2(planck_x,planck)
int3 = integral3(planck_x,planck)

print
print "Stephan-Boltzman:", total
print "First method:",int1*const.pi, ", diffence:", abs(total-int1*const.pi), abs(total-int1*const.pi)/total
print "Second method:",int2*const.pi, ", diffence:", abs(total-int2*const.pi), abs(total-int2*const.pi)/total
print "Third method:",int3*const.pi, ", diffence:", abs(total-int3*const.pi), abs(total-int3*const.pi)/total
print



matplotlib.rcParams.update({'font.size': 12, 'font.family': 'serif'})



figure(figsize=(8,5))
title('Planck')
xlabel('$\lambda$ (nm)')
ylabel('B$_{\lambda}$(T) (W.m$^{-2}$.nm$^{-1}$.sr$^{-1}$)')
xlim([0,5000])
grid()
plot(planck_x,planck)
savefig("planck.png")
