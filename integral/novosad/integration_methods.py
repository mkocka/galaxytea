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
	  