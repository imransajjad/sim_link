#!/usr/bin/python

"""Some ode solvers.

Maintained on Github
https://github.com/imransajjad/sim_link

by
Imran Sajjad
imransajjad89@gmail.com

t is considered the independent variable or
integration wrt t

For all solver functions, the order of input arguments should be as follows:

f, arg2 to f, arg3 to f, ..., argN to f, Trange, init X

where dx/dt = f(t,x,arg2,arg3,...argN)

arg0 and arg1 are t and x respectively and will be handled by the funtion

prefrably init X and Trange should be np.arrays
and f should also return an np.array the size of init X

Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. 
"Runge-Kutta Method" and "Adaptive Step Size Control for Runge-Kutta." 
16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific 
Computing, 2nd ed. Cambridge, England: Cambridge University Press, 
pp. 704-716, 1992.

"""
import time
from math import sqrt



def integrate(f,*args):
	""" simple dirty integration
	fixed step:	x = x+f*dt
	it is expected that f(t,x,*args)
	the order of input arguments should be as follows #
	function, arg1 to f, arg2 to f, ..., argn to f, Trange, init X"""

	x = args[-1]
	X = [x]
	T = [args[-2][0]]

	for t in args[-2]:

		k = f(t,x, *args[0:-2] )
		x = x + k*dt
		X.append(x)
		T.append(t)
	return T,X

def trapz(f,*args):
	""" trapezoidal integration
	fixed step:	x = x+0.5*(f(x)+f(x+f(x)*dt))*dt
	it is expected that f(t,x,*args)
	the order of input arguments should be as follows #
	function, arg1 to f, arg2 to f, ..., argn to f, Trange, init X"""

	x = args[-1]
	X = [x]
	T = [args[-2][0]]

	dt = args[-2][1]-args[-2][0]
	for t in args[-2]:
		k1 = f(t, x, *args[0:-2] )
		k2 = f(t+dt, x+k1*dt, *args[0:-2] )
		x = x +  0.5*(k1+k2)* dt

		X.append(x)
		T.append(t)
	return T,X



def rungekutta4(f, *args):
	""" trapezoidal integration
	fixed step:	four terms of rungekutta4
	it is expected that f(t,x,*args)
	the order of input arguments should be as follows #
	function, arg1 to f, arg2 to f, ..., argn to f, Trange, init X"""

	x = args[-1]
	X = [x]
	T = [args[-2][0]]

	dt = args[-2][1]-args[-2][0]
	for t in args[-2]:
		
		k1 = f(t, x, *args[0:-2])
		k2 = f(t+dt/2, x+k1*dt/2, *args[0:-2])
		k3 = f(t+dt/2, x+k2*dt/2, *args[0:-2])
		k4 = f(t+dt, x+k3*dt, *args[0:-2])
		
		x = x+ (k1+2*k2+2*k3+k4)/6*dt
		X.append(x)
		T.append(t+dt)
	return T,X

def rungekutta4ad(f, *args, **kwargs):
	
	t = args[-2][0]
	dt = args[-2][1]-args[-2][0]
	Tf = args[-2][-1]
	x = args[-1]

	T = [t]
	X = [x]

	P = {"min_dt": 1e-5, "e_tol": 1e-6, "gain": 0.04, "adaptive": False,\
	 "realtime": False, "outcall": None, "plotcalls": [], "plottime" : 0.001}
	P.update(kwargs)

	O = 7
	C = [0.0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1.0, 1.0]

	A = [ [],\
	[1.0/5],\
	[3.0/40, 9.0/40],\
	[44.0/45, -56.0/15, 32.0/9],\
	[19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729],\
	[9017.0/3168, -355.0/33, -46732.0/5247, 49.0/176, -5103.0/18656],\
	[35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84] ]

	Bup = [35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0.0]
	Bdn = [5179.0/57600, 0.0, 7571.0/16695, 393.0/640, -92097.0/339200, 187.0/2100, 1.0/40]

	k = [0.0]*O

	# O = 4
	# C = [0.0, 0.5, 0.5, 1.0]
	# A = [ [], [0.5], [0.0, 0.5], [0.0, 0.0, 1.0]]
	# Bup = [1.0/6, 1.0/3, 1.0/3, 1.0/6]
	# Bdn = [1.0/6, 1.0/3, 1.0/3, 1.0/6]
	# k = [0.0]*O

	if not P["realtime"]:
		Ttarget = Tf
	else:
		sys_time = time.time()
		plot_time = time.time()
		ti = 1
		len_T = len(args[-2])
		Ttarget = args[-2][ti]

	Y = []
	if P["outcall"]:
		Y.append( P["outcall"](t,x, *args[0:-2]) )

	print("ode_solver with config")
	for p in P.keys():	print(p+":", P[p])
	print("")

	while t < Ttarget:
		
		for i,(c,a_row) in enumerate(zip(C,A)):
			k[i] = f(t+c*dt, x + dt*sum([k[j]*a for j,a in enumerate(a_row)])  , *args[0:-2])
		
		gradup = sum([k[i]*b for i,b in enumerate(Bup)])

		if P["adaptive"]:
			graddn = sum([k[i]*b for i,b in enumerate(Bdn)])
			E =gradup-graddn
			E = sqrt(sum([e*e for e in E]))

			margin = float(E/P["e_tol"])
			# print margin, t
			
			if dt < P["min_dt"] or margin < 1.0:
				# append result if error within bounds
				t = t+dt
				x = x+ gradup*dt
				T.append(t)
				X.append(x)
				dt = (1.0+P["gain"])*dt
				if P["outcall"]:
					Y.append( P["outcall"](t,x, *args[0:-2]) )
			else:
				# retry
				dt = (1.0-P["gain"])*dt
				# print "retrying with dt = ", dt, margin
		else:
			t = t+dt
			x = x+ gradup*dt
			T.append(t)
			X.append(x)
			if P["outcall"]:
				Y.append( P["outcall"](t,x, *args[0:-2]) )

		if P["realtime"]:
			ti +=1
			if ti < len_T:
				Ttarget = args[-2][ti]
				cur_time = time.time()
				# print "here plotting"
				# print (args[-2][ti]-args[-2][ti-1]) - (time.time()-sys_time) > 0.001
				if (args[-2][ti]-args[-2][ti-1]) - (time.time()-sys_time) > P["plottime"]:
					plot_time = time.time()
					if P["plotcalls"]:
						# print "here2"
						for p in P["plotcalls"]:
							p(T,X,* ((Y,) if P["outcall"] else ()) )
					# print time.time()-plot_time
						
						# cur_time = time.time()
				# print time.time()-sys_time , (args[-2][ti]-args[-2][ti-1])
				while (time.time()-sys_time) < (args[-2][ti]-args[-2][ti-1]):
					# print "waiting"
					pass

					# cur_time = time.time()
				sys_time = time.time()
			
			

			
		

	if P["plotcalls"] and t >= Tf:
		for p in P["plotcalls"]:
			p(T,X,* ( (Y,) if P["outcall"] else () ) )

	if not P["outcall"]:
		return T,X
	else:
		return T,X,Y



def test_f(t,x):
	dx = np.array([0.0]*2)

	dx[0] = x[1]
	dx[1] = -2*x[0] - 0.5*x[1]

	return (dx,0.0)

def normpdf(x,mu,sigma):
	xad = (x-mu)/sigma
	return np.exp(-xad**2/2)/np.sqrt(2*np.pi)/sigma

def fX(x):
	mu = 1.0
	sigma = 0.1
	return normpdf(x,mu,sigma)

def fY(x,X):
	mu = X
	sigma = 0.01
	return normpdf(x,mu,sigma)

def test1():
	T = np.arange(0,10,0.05)
	x0 = np.array([2.0,0])

	k = 2
	m = 1
	b = 0.5
	wn = np.sqrt(k/m)
	gamma = b/2/np.sqrt(k*m)
	wd = wn*np.sqrt(1-gamma**2)
	sigma = gamma*wn
	s1 = np.array([-sigma + wd*1.0j])
	s2 = np.array([-sigma - wd*1.0j])
	print(k,m,b, gamma, wd, sigma, wn)

	x00 = 2
	v00 = 0
	alpha = x00/2
	beta = (v00+sigma*x00)/2/wd
	
	X0 = [ 1*np.exp(s1*t) + 1*np.exp(s2*t) for t in T]
	X0 = [ 2*(alpha*np.cos(wd*t) + \
		beta*np.sin(wd*t)) *np.exp(-sigma*t) for t in T]
	X1 = integrate(test_f, T, x0)
	X2 = trapz(test_f, T, x0)
	X3 = rungekutta4(test_f, T, x0)

	plt.plot(T, X0, label='true')
	plt.plot(T, X1[1][:,0], label='int')
	plt.plot(T, X2[1][:,0], label='trapz')
	plt.plot(T, X3[1][:,0], label='rk4')
	plt.legend()
	plt.show()

def test2():
	x = np.arange(-50,50,0.01)

	e = lambda x,u : np.exp(-(x/u)**2)/u
	D = lambda x,y,u1,u2 : \
		(np.array([np.sqrt(2/np.pi*u1*u2)*e(x,u1)*e(x,u2)]),0.0)

	u1 =10.0
	U2 = np.arange(0.1,2.51,0.05)
	Y = 0*U2
	for i,u2 in enumerate(U2):
		print(u2)
		_,y,_ = rungekutta4( D,  u1, u2, x, np.array([0.0]))
		Y[i] = y[-1][0]

	plt.plot(U2,Y)
	plt.plot(U2,np.sqrt(2/(U2/u1 + u1/U2)))
	plt.legend()
	plt.show()


if __name__ == '__main__':
	import numpy as np
	import matplotlib.pyplot as plt

	print(__doc__)
	test1()