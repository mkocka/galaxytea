
import matplotlib.pyplot as plt


def galaxy_tea_1(n):
	Rhv = 1
	alpha = 0.4
	M16 = 5
	m1 = 12
	Sigma = []
	H = []
	rho = []
	Tc = []
	tau = []
	ny = []
	vr = []
	R_10 = []
	for R in range(Rhv+1,n,1):
		R_10.append(float(R))	
		R10 = (float(R))
		f = (1-(Rhv/float(R))**0.5)**0.25
		
		Sigma.append(5.2*alpha**(-4.0/5)*M16**(7.0/10)*m1**(1/4)*R10**(-3/4)*f**(14.0/5))
		H.append(1.7*10**(8)*alpha**(-1.0/10)*M16**(3.0/20)*m1**(-3.0/8)*R10**(9.0/8)*f**(3.0/5))
		rho.append(3.1*10**(-8)*alpha**(-7.0/10)*M16**(11.0/20)*m1**(5/8)*R10**(-15/8)*f**(11.0/5))
		Tc.append(1.4*10**4*alpha**(-1.0/5)*M16**(3.0/10)*m1**(1.0/4)*R10**(-3/4)*f**(6.0/5))
		tau.append(190*alpha**(-4.0/5)*M16**(1.0/5)*f**(4.0/5))
		ny.append(1.8*10**(14)*alpha**(4.0/5)*M16**(3.0/10)*m1**(-1/4)*R10**(3.0/4)*f**(6.0/5))
		vr.append(2.7*10**(4)*alpha**(4.0/5)*M16**(3.0/10)*m1**(-1/4)*R**(-1.0/4)*f**(-14.0/5))

		
	fig = plt.figure()
	plt.title("Sigma na R_10")
	plt.xlabel('R_10')
	plt.ylabel('Sigma')
	plt.plot(R_10,Sigma,'r')
	plt.savefig("Sigma_0.4.jpg")
	plt.close(fig)
	
	fig = plt.figure()
	plt.title("H na R_10")
	plt.xlabel('R_10')
	plt.ylabel('H')
	plt.plot(R_10,H,'r')
	plt.savefig("H_0.4.jpg")
	plt.close(fig)
	
	fig = plt.figure()
	plt.title("rho na R_10")
	plt.xlabel('R_10')
	plt.ylabel('rho')
	plt.plot(R_10,rho,'r')
	plt.savefig("rho_0.4.jpg")
	plt.close(fig)
	
	fig = plt.figure()
	plt.title("Tc na R_10")
	plt.xlabel('R_10')
	plt.ylabel('Tc')
	plt.plot(R_10,Tc,'r')
	plt.savefig("Tc_0.4.jpg")
	plt.close(fig)
	
	fig = plt.figure()
	plt.title("tau na R_10")
	plt.xlabel('R_10')
	plt.ylabel('tau')
	plt.plot(R_10,tau,'r')
	plt.savefig("tau_0.4.jpg")
	plt.close(fig)
	
	fig = plt.figure()
	plt.title("ny na R_10")
	plt.xlabel('R_10')
	plt.ylabel('ny')
	plt.plot(R_10,ny,'r')
	plt.savefig("ny_0.4.jpg")
	plt.close(fig)
	
	fig = plt.figure()
	plt.title("vr na R_10")
	plt.xlabel('R_10')
	plt.ylabel('vr')
	plt.plot(R_10,vr,'r')
	plt.savefig("vr_0.4.jpg")
	plt.close(fig)

galaxy_tea_1(100)	
		
