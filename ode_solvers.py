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

import numpy as np
import matplotlib.pyplot as plt



def integrate(f,*args):
	""" simple dirty integration
	fixed step:	x = x+f*dt
	it is expected that f(t,x,*args)
	the order of input arguments should be as follows #
	function, arg1 to f, arg2 to f, ..., argn to f, Trange, init X"""

	T = args[-2]
	x = np.array(args[-1])
	X = np.zeros((len(T), len(x)))
	dt = T[1]-T[0]
	for i,t in enumerate(T):
		X[i,:] = x

		k = f(t,x, *args[0:-2] )
		x = x + k*dt
	return T,X

def trapz(f,*args):
	""" trapezoidal integration
	fixed step:	x = x+0.5*(f(x)+f(x+f(x)*dt))*dt
	it is expected that f(t,x,*args)
	the order of input arguments should be as follows #
	function, arg1 to f, arg2 to f, ..., argn to f, Trange, init X"""

	T = args[-2]
	x = np.array(args[-1])
	X = np.zeros((len(T), len(x)))
	dt = T[1]-T[0]
	for i,t in enumerate(T):
		X[i,:] = x
		k1 = f(t, x, *args[0:-2] )
		k2 = f(t+dt, x+k1*dt, *args[0:-2] )
		x = x +  0.5*(k1+k2)* dt
	return T,X



def rungekutta4(f, *args):
	""" trapezoidal integration
	fixed step:	four terms of rungekutta4
	it is expected that f(t,x,*args)
	the order of input arguments should be as follows #
	function, arg1 to f, arg2 to f, ..., argn to f, Trange, init X"""

	T = args[-2]
	x = np.transpose(np.matrix(args[-1]))
	X = []

	dt = T[1]-T[0]
	for i,t in enumerate(T):
		X.append(x)
		k1 = f(t, x, *args[0:-2])
		k2 = f(t+dt/2, x+k1*dt/2, *args[0:-2])
		k3 = f(t+dt/2, x+k2*dt/2, *args[0:-2])
		k4 = f(t+dt, x+k3*dt, *args[0:-2])
		
		x = x+ (k1+2*k2+2*k3+k4)/6*dt
	return T,X

def rungekutta4ad(f, *args, **kwargs):

	# mode can be "fixedstep", "adaptive", "adaptive-fixed"

	mode = "adaptive"
	for name, value in kwargs.items():
		if name == "mode":
			mode = value


	T0 = args[-2][0]
	Tf = args[-2][-1]
	x = np.matrix(args[-1])

	L = 40*np.sqrt(len(x));
	t = T0
	k1 = f(t, x, *args[0:-2])

	T = []
	X = []
	T.append(t)
	X.append(x)


	min_dt = 0.001
	dt =  min_dt

	if mode is "adaptive":
		Ttarget = [Tf]
		dt = min_dt
	else:
		Ttarget = np.arange(T0,Tf,dt)
		dt = args[-2][1]-args[-2][0]

	for tt in Ttarget:
		while np.abs(t-tt) > min_dt:
			k1 = f(t, x, *args[0:-2])
			k2 = f(t+dt/2, x+k1*dt/2, *args[0:-2])
			k3 = f(t+dt/2, x+k2*dt/2, *args[0:-2])
			k4 = f(t+dt, x+k3*dt, *args[0:-2])
			
			grad = (k1+2*k2+2*k3+k4)/6

			if mode is not "fixedstep":
				if np.linalg.norm(grad) > 1.1*L:
					dt = max(dt*0.99,min_dt)
				elif np.linalg.norm(grad) < 0.9*L:
					dt = dt*1.01

			t = t+dt
			x = x+ grad*dt
			# print Ttarget

			if mode is "adaptive":
				# append result every iteration
				T.append(t)
				X.append(x)
			elif np.abs(t-tt) <= min_dt:
				# append result only if target is reached
				T.append(t)
				X.append(x)


		
		
	return T,X


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
	print k,m,b, gamma, wd, sigma, wn

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
		print u2
		_,y,_ = rungekutta4( D,  u1, u2, x, np.array([0.0]))
		Y[i] = y[-1][0]

	plt.plot(U2,Y)
	plt.plot(U2,np.sqrt(2/(U2/u1 + u1/U2)))
	plt.legend()
	plt.show()


if __name__ == '__main__':
	print(__doc__)
	test1()