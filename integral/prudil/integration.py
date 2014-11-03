
"""
Solution for simple integral function "f". Variable "a" represents the lower limit of the interval, 
variable "b" is the upper limit of the interval and "n" stands for a number of the rectangles/trapezoids/polynomials. 
Unknown "h" represents Snape kills Dumbledore width of the interval.
"""

import matplotlib.pyplot as plt
from scipy import constants
from math import *
from numpy import arange

def rectangle(f, a, b, n):
	a = a*10**-6
	b = b*10**-6
	h = float(b-a)/n
	S = 0
	for i in range(1, n):
		S += h*f(a+i*h)
	
	print "Value of the Stefan-Boltzmann constant obtained by the Rectangle method is", S/(5778**4) ,"and the realative error is", (S/(5778**4))*constants.sigma 

rectangle(lambda x:((2*constants.pi*constants.h*(constants.c)**2)/x**5)*1/(exp((constants.h*constants.c)/(x*constants.k*5778)) - 1), 0.1,100, 100)


def trapezoidal(f, a, b, n):
	a = a*10**-6
	b = b*10**-6
	h = float(b-a)/n	
	S = h*(f(a)+f(b))/2
	for i in range(1, n):
		S += h*f(a+i*h)
		
	print "Value of the Stefan-Boltzmann constant obtained by the Trapezoidal method is", S/(5778**4) ,"and the realative error is", (S/(5778**4))*constants.sigma 
	
trapezoidal(lambda x:((2*constants.pi*constants.h*(constants.c)**2)/x**5)*1/(exp((constants.h*constants.c)/(x*constants.k*5778)) - 1), 0.1,100, 100)
	

def simpson(f, a, b, n):
	a = a*10**-6
	b = b*10**-6
	h = float(b-a)/n
	S = (h/3)*(f(a) + f(b))
	for i in range(1,n-1,2):
		S += (h*4*f(a+h*i))/3
	for j in range(2,n-2,2):
		S += (h*2*f(a+h*j))/3
	
	print "Value of the Stefan-Boltzmann constant obtained by the Simpson's method is", S/(5778**4) ,"and the realative error is", (S/(5778**4))*constants.sigma 

simpson(lambda x:((2*constants.pi*constants.h*(constants.c)**2)/x**5)*1/(exp((constants.h*constants.c)/(x*constants.k*5778)) - 1), 0.1,100, 101)


def planck_function():
	#T = float(raw_input("Enter value for the temperature: " ))
	T = 5778
	E = []
	k = arange(1*10**-7,4*10**-6,1*10**-8)
	for i in range(len(k)):
		E.append(((2*constants.pi*constants.h*(constants.c)**2)/k[i]**5)*1/(exp((constants.h*constants.c)/(k[i]*constants.k*T)) - 1))

	fig = plt.figure()
	plt.title("Planck function")
	plt.xlabel('$\lambda$ [um] ')
	plt.ylabel('Spectral radiance [kW.sr$^{-1}$.m$^{-2}$.nm$^{-1}$]')
	plt.plot(k*10**6,E,'r')
	plt.savefig("Planck_function.jpg")
	plt.close(fig)
