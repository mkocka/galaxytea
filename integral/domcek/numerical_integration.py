#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt 
import matplotlib
from scipy import constants as const
from pylab import *

#

def planck(T,lambda_i,lambda_f,step=10**(-1)):
	x ,y = ([] for i in range(2))
	for i in arange(lambda_i,lambda_f,step): # i*10**(-9) - correctuin to nm
		f = 10**(-9)*(2*const.h*const.c**2)/(((i*10**(-9))**5) *(exp((const.c*const.h)/((i*10**(-9))*const.k*T))-1))
		y.append(f)
		x.append(i)
	return x,y
	
def plot_f(x,y): #plots one planck function
	plt.plot(x,y)
	plt.xlabel('$\lambda$ [nm]')
	plt.ylabel('Spectral radiance [W.sr$^{-1}$.m$^{-}$.nm$^{-1}$')
	plt.grid()	
	plt.savefig("planck")
	plt.gcf().clear()	
	
#first 3 parameters characterize range of temperatures and step between them
#other 2 parameters characterize range of planck function
def plot_multiple(T_i,T_f,step,lambda_i,lambda_f): #plots multiple planck into one graph
	for i in range(T_i,T_f,step):
		a, b = planck(i,lambda_i,lambda_f)
		plt.plot(a,b)
	plt.xlabel('$\lambda$ [nm]')
	plt.ylabel('Spectral radiance [W.sr$^{-1}$.m$^{-}$.nm$^{-1}$')
	plt.grid()	
	plt.savefig("planck_functions")
	plt.gcf().clear()

def midpoint_rule(a,b): #rectangle rule
	dx = a[1]-a[0]
	summation=0
	for i in range(1,len(a)): #first and last value excluded
		summation += dx*b[i]
	return summation
	
def trapezoid_rule(a,b):
	dx = a[1]-a[0]
	summation=0.5*b[0]+0.5*b[len(b)-1]
	for i in range(1,len(a)): #first and last value excluded
		summation += dx*b[i]
	return summation	

def simpson(a,b):
	dx = a[1]-a[0]
	summation = dx*(b[0]+b[len(b)-1])/3
	ax ,bx = ([] for i in range(2))
	if (len(a)%2==0): # this rule needs odd number of values
		a.pop() and b.pop()
	for i in range(len(a)):
		if i%2==0: summation +=4*b[i]*dx/3
		if i%2==1: summation +=2*b[i]*dx/3
	return summation


plot_multiple(5000,10000,1000,10,2000)


#defining boundries of calculating function
lambda_initial=10
lambda_final=3000

#temperature of blackbody
temperature=5778

#calculating planck function
a, b = planck(temperature,lambda_initial,lambda_final) 
plot_f(a,b)
 
#theoretical stefan - boltzman value
stef_boltz = const.sigma*temperature**4


print 30*"-"
print "Planck function parameters"
print 30*"-"
print "Temperature:", temperature
print "Range of integration in nm: (", lambda_initial, ",", lambda_final,")"
print
print "Calculated values:"
print 30*"-"
#integration gives us value dependent on sr**-1, need to mutliply by appropriate constant (pi)
#for more information: http://en.wikipedia.org/wiki/Planck's_law#Stefan.E2.80.93Boltzmann_law
print "Stephan-Boltzman law:", stef_boltz 
print "Midpoint rule method:", midpoint_rule(a,b)*const.pi
print "Trapezoid rule method:", trapezoid_rule(a,b)*const.pi
print "Simpson's rule method:", simpson(a,b)*const.pi
print
print "Relative errors"
print 30*"-"
print "Midpoint rule method:", (stef_boltz  -  midpoint_rule(a,b)*const.pi) / stef_boltz
print "Trapezoid rule method:",  (stef_boltz  - trapezoid_rule(a,b)*const.pi) / stef_boltz
print "Simpson's rule method:", (stef_boltz  - simpson(a,b)*const.pi) / stef_boltz
