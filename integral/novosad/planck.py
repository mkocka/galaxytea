import matplotlib.pyplot as plt 
import matplotlib
import numpy as np
import math
import scipy
from integration_methods import int_rectangular  
from integration_methods import int_trapezoid
from integration_methods import int_simpson

##### T ... temperature [K]
##### a ... left point of interval ######## a and b are wave lenghts
##### b ... right point of interval
##### N ... length of step


###############################    planck law #########################################################
def planck(T,a,b,N):
      planck = []
      h = 6.636 
      c = 3.0 
      k = 1.38 
      for i in np.arange(a,b,N):
	    exp = np.exp(h * c * 10**(-3)/(i * 10**(-9) * k * T))
	    planck.append((2.0 * h * 10**(-34) * c**2 * 10**(16) * 10**(45))/(i**5 * (exp-1)))
      return planck
  
###############################  imput to planck law ############################################
a = 10.0   ### in nm
b = 3000.0
N = 1000
h = (b-a)/N


x = []
for i in np.arange(a,b,h):
      x.append(i)

###############################  plot        ###################################################### 
fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(x,planck(4000.0,a,b,h),label=r"4000 K")
axes.plot(x,planck(4500.0,a,b,h),label=r"4500 K")
axes.plot(x,planck(5000.0,a,b,h),label=r"5000 K")
axes.plot(x,planck(5500.0,a,b,h),label=r"5500 K")
axes.plot(x,planck(5778.0,a,b,h),label=r"5778 K")
axes.plot(x,planck(6000.0,a,b,h),label=r"6000 K")
axes.set_xlabel(r'$\lambda \, [\mathrm{nm}]$')
axes.set_ylabel(r'$B_{\lambda} \, [\mathrm{W \cdot sr^{-1} \cdot m^{-3} }]$')
axes.set_title('Planck law')
axes.legend(loc=1)
fig.savefig("planck.png")

###########################   integration  ##########################################################

sol1=int_rectangular(planck(5778.0,a,b,h),a,b,N)*10**(-9)*3.14
sol2=int_trapezoid(planck(5778.0,a,b,h),a,b,N)*10**(-9)*3.14
sol3=int_simpson(planck(5778.0,a,b,h),a,b,N)*10**(-9)*3.14
print(sol1," rectangular")
print(sol2," trapezoid")
print(sol3," Simpson")

sol_analytical = 5.670373 * 10**(-8) * 5778.0**4
print(sol_analytical," Stefan-Boltzmann law")
