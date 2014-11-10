import matplotlib.pyplot as plt 
import matplotlib
import numpy as np
import math
import scipy


### fct - function
### a   - left point of interval
### b   - right point of interval
### N   - number of steps ... for Simpson rule N must be odd!!!!
 
#####################    declaration of integration methods  #####################

def int_rectangular(fct,a,b,N):
	  h=(b-a)/N
	  rectangular = 0.0
	  for i in range(N):
		rectangular = rectangular + fct[i]*h
	  return rectangular
	
def int_trapezoid(fct,a,b,N):
	  h=(b-a)/N
	  trapezoid = (fct[0] + fct[N-1])*h/2
	  for i in range(1,N-1):
		trapezoid = trapezoid + fct[i]*h
	  return trapezoid
	
def int_simpson(fct,a,b,N):
	  h=(b-a)/N
	  simpson = (fct[0] + fct[N-1])*h/3
	  for i in np.arange(1,N-1,2):
		simpson = simpson + 4*h*fct[i]/3
	  for i in np.arange(2,N-2,2):
		simpson = simpson + 2*h*fct[i]/3
	  return simpson	
	  

################################################################################################
################################## Parameters of compact object ################################
alpha = 0.5			#parameter of accretion [something]
M = 1.0				#change of mass of compact object [[10**16 g * s**(-1)]]
m = 5.0				#mass of compact object [M_sun]
R_star = 10.0**(-4) 		#radius of compact object [10**10 cm = 100 000 km, so 10**(-5)= 1 km]		10 km
R_out = 10.0**(-2)		#outer radius of disc [10**10 cm = 100 000 km, so 10**(-5)= 1 km]		1000 km
R_step = 10.0**(-6) 		#steps in computation [10**10 cm = 100 000 km, so 10**(-5)= 1 km]		100 m 

##################################################################################################
T_c = []
Rx = []


for R in np.arange(R_star+R_step,R_out,R_step):
	Rx.append((R_star + R_step +R) * 10**5) 									# km
	f=(1-(R_star/R)**(1.0/2))**(1.0/4)										# no dimension
	T_c.append(1.4 * 10**4 * alpha**(-1.0/5) * M**(3.0/10) * m**(1.0/4) * R**(-3.0/4) * f**(11.0/5))		# K


############################### spectrum   #####################################################

##############################################################################################################
##############################################################################################################
r_in = 1.0
theta = 1.0
D = 1.0
##############################################################################################################
##############################################################################################################
T_in = max(T_c)
size = len(T_c)    
T_out = T_c[size-1]
index_of_T_in = T_c.index(max(T_c))


h = 6.636 
c = 3.0 
k = 1.38 
energy = []
new_T_c = []
integral = []
for i in range(index_of_T_in,size):
	new_T_c.append(T_c[i])
		    
		    
for E in np.arange(20.0,2000.0,10):								#### decleration of interval of energies
	energy.append(E)
	planck_for_constant_energy = []
	fct = []
	
	for i in range(len(new_T_c)):
	      planck_for_constant_energy.append(2.0 * 1.6**3 * E**3 * 10**(29)/( h**3 * c**2  * (np.exp(1.6 * E/(k * T_c[i] * 10**(-4)))-1))) 
	      fct.append((new_T_c[i]/(T_in))**(-11.0/3)*planck_for_constant_energy[i] /(T_in))
		    
	integral.append(int_rectangular(fct,T_out,T_in,len(new_T_c)))
f_d = []
for i in range(len(integral)):
	f_d.append(8 * 3.14 * math.cos(theta) * r_in**2 * integral[i]/( 3 * D**2))

###############################  plots ##########################################################
fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,T_c)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$T_c \, [\mathrm{K}]$')
axes.set_title('Temperature')
fig.savefig("T_c.png")

### plot of stpectrum

fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(energy,f_d)
axes.set_xlabel('$E \, [\mathrm{eV}]$')
axes.set_ylabel('$f_d \, [\mathrm{neco}]$')
axes.set_title('Spectrum')
fig.savefig("Full_Spectrum.png")


### plot of fct to integration
#fig = plt.figure()
#axes = fig.add_axes([0.17,0.17,0.75,0.75])
#axes.plot(new_T_c,fct,label=r"max K")
#axes.set_xlabel(r'$T \, [\mathrm{K}]$')
#axes.set_ylabel('neco')
#axes.set_title('Spectrum for one energy')
#axes.legend(loc=1)
#fig.savefig("spectrum_for_one_energy.png")