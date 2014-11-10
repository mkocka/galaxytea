#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt 
import matplotlib
from scipy import constants as const
from pylab import *
from numpy import *


def solutions(r_in,r_out,step,alfa,M_16,m_1,R_hv):
	#defining lists
	list_function = arange(r_in,r_out,step)
	R_10_l, temperature_l= ([] for i in range(2))
	#computation and appending to lists
	for R_10 in list_function:
		f=(1-((R_hv)/(R_10))**(1.0/2))**(1.0/4)	
		temperature = 1.4*10**4*alfa**(-1.0/5)*M_16**(3.0/10)*m_1**(1.0/4)*R_10**(-3.0/4)*f**(6.0/5)
		temperature_l.append(temperature)
		R_10_l.append(R_10)
	return R_10_l, temperature_l
	
def midpoint_rule(a,b): #rectangle rule
	dx = a[0]-a[1] #changed sign to - (temperature is decreasing)
	summation=0
	for i in range(1,len(a)): #first and last value excluded
		summation += dx*b[i]
	return summation
	
def model(T_in,T_out,temperaturestep,E_initial,E_final,energystep,R_in,theta,distance):
	temperature = arange(T_in,T_out,-temperaturestep)
	energy = arange(E_initial,E_final,energystep)
	C = ((8*const.pi*R_in*cos(radians(theta)))/(3*distance**2))
	x, y ,logx = ([] for i in range(3))
	for j in energy:
		T ,f = ([] for i in range(2))
		for i in temperature:
			planck = (2*(j*c_e)**3*10**29)/((c_h**3*c_c**2) *(exp((j*const.e)/(const.k*i))-1)) 
			#energy in planck multipied by eV, h**3 in code = 0, instead using 10**29 gain from e**3/h**3.c**2
			T.append(i)
			f.append((float(i)/T_in)**(-11.0/3)*planck/T_in)
		x.append(j)
		logx.append(log10(j))
		y.append(C*midpoint_rule(T,f))
		print j
	return x,y,logx
	
#parameters for temperature profile of the disk 
r_in =1.0001*10**(-4) #r_in is radius of lower boundry of integration
r_out =10**(-2)
step = 10**(-6)
alfa = 0.5
M_16 = 63
m_1 = 1.5
R_hv = 1.0*10**(-4)
	
R,T = solutions(r_in,r_out,step,alfa,M_16,m_1,R_hv)

indexmax = argmax(T) #index of maximum temperature in list 
T_in= T[indexmax] 
R_in= R[indexmax] #R_in is radius of maximum temperature

#parameters for function model	
c_h = const.h*10**(34)
c_c = const.c*10**(-8)
c_e = const.e*10**19
parsec = 3.1*10**18 #in centimeters

E_initial=1
E_final=20000
energystep = 100
T_out = 10**3
temperaturestep = 10**3
theta = 0
distance= 8*parsec
#

x,y,logx = model(T_in,T_out,temperaturestep,E_initial,E_final,energystep,R_in,theta,distance)

			
plt.plot(x, y)
plt.title('Blackbody Spectrum')
plt.xlabel('E [eV]')
plt.ylabel('Intensity ')
plt.grid()
plt.savefig("MCD")
plt.gcf().clear()		

plt.plot(logx, y)
plt.title('Blackbody Spectrum')
plt.xlabel('log E [eV]')
plt.ylabel('Intensity W.m$^{-2}$.eV ')
plt.grid()
plt.savefig("MCDlog")
plt.gcf().clear()	
