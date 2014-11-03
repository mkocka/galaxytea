import matplotlib.pyplot as plt 
import matplotlib
import numpy as np
import math

alpha = 0.5			#parameter of accretion [something]
M = 1.0				#change of mass of compact object [[10**16 g * s**(-1)]]
m = 5.0				#mass of compact object [M_sun]
R_star = 10.0**(-4) 		#radius of compact object [10**10 cm = 100 000 km, so 10**(-5)= 1 km]		10 km
R_out = 10.0**(-2)		#outer radius of disc [10**10 cm = 100 000 km, so 10**(-5)= 1 km]		1000 km
R_step = 10.0**(-6) 		#steps in computation [10**10 cm = 100 000 km, so 10**(-5)= 1 km]		100 m 



Sigma = []
H = []
rho = []
T_c = []
tau = []
nu = []
v_R = []
Rx = []

for R in np.arange(R_star+R_step,R_out,R_step):
	Rx.append((R_star + R_step +R) * 10**5) 									# km
	f=(1-(R_star/R)**(1.0/2))**(1.0/4)										# no dimension
	Sigma.append(5.2 * alpha**(-4.0/5) * M**(7.0/10) * R**(-3.0/4) * f**(14.0/5))					# g cm**-2
	H.append(1.7 * 10**8 * alpha**(-1.0/10) * M**(3.0/20) * m**(-3.0/8) * R**(9.0/8) * f**(3.0/5))			# cm
	rho.append(3.1 * 10**(-8.0) * alpha**(-7.0/10) * M**(11.0/20) * m**(5.0/8) * R**(-15.0/8) * f**(11.0/5))	# g cm**-3
	T_c.append(1.4 * 10**4 * alpha**(-1.0/5) * M**(3.0/10) * m**(1.0/4) * R**(-3.0/4) * f**(11.0/5))		# K
	tau.append(190 * alpha**(-4.0/5) * M**(1.0/5) * f**(4.0/5))							# no dimension
	nu.append(1.8 * 10**14 * alpha**(4.0/5) * M**(3.0/10) * m**(-1.0/4) * R**(3.0/4) * f**(6.0/5))			# cm**2 s**-1
	v_R.append(2.7 * 10**4 * alpha**(4.0/5) * M**(3.0/10) * m**(-1.0/4) * R**(-1.0/4) * f**(-14.0/5))		# cm s**-1
	
fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,Sigma)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$\Sigma \, [\mathrm{g \cdot cm^{-2} }]$')
axes.set_title('Surface density')
fig.savefig("Sigma.png")

fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,H)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$H \, [\mathrm{cm}]$')
axes.set_title('Height')
fig.savefig("H.png")

fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,rho)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$ \\rho \, [\mathrm{g \cdot cm^{-3}] }$')
axes.set_title('Density')
fig.savefig("density.png")

fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,T_c)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$T_c \, [\mathrm{K}]$')
axes.set_title('Temperature')
fig.savefig("T_c.png")

fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,tau)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$\\tau$')
axes.set_title('Opacity')
fig.savefig("tau.png")

fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,nu)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$\\nu \, [\mathrm{cm \cdot s^{-1} }]$')
axes.set_title('Viscosity')
fig.savefig("nu.png")

fig = plt.figure()
axes = fig.add_axes([0.17,0.17,0.75,0.75])
axes.plot(Rx,v_R)
axes.set_xlabel('$R \, [\mathrm{km}]$')
axes.set_ylabel('$v_R \, [\mathrm{cm \cdot s^{-1} }]$')
axes.set_title('Radial speed')
fig.savefig("v_R.png")


