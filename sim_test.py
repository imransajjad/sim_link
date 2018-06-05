import numpy as np
import ode_solvers as ode
import sim_link as sl
import sim_lib as sll

import matplotlib
# matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt

def B():

	def der(t,x,b):
		xdot = np.matrix([])
		return xdot
	def out(t,x,b):
		y = 7*b[0,0]
		return y

	B.der = der
	B.out = out
	B.namestring = "function B"

	B.inargs = 1
	B.cstates = 0
	B.passargs = [0]

def A():

	def der(t,x,a):
		xdot = np.matrix([-4.0*x[0,0]])
		return xdot
	def out(t,x,a):
		y = np.matrix([x[0,0]+a[0,0]])
		return y

	A.der = der
	A.out = out
	A.namestring = "function A"

	A.inargs = 1
	A.cstates = 1
	A.passargs = [0]

def D():

	def der(t,x,d1,d2):
		xdot = -x
		return xdot
	def out(t,x,d1,d2):
		# print t, x, d1, d2
		y = np.matrix([x[0,0]+d1[0,0]+d2[0,0]])
		return y

	D.der = der
	D.out = out
	D.namestring = "function D"

	D.cstates = 2
	D.inargs = 2
	D.passargs = [0,1]





def G():

	def der(t,x,u):
		# x1dot = x2dot
		# x2dot = -1.0*x1 -2.0*x2 + u
		# xdot = np.array([ x[1], -1.0*x[0] -2.0*x[1] + u[0] ])
		A = np.matrix('0 1; -1 -2')
		B = np.matrix('0;1')
		xdot = A*x + B*u
		return xdot
	def out(t,x,u):
		# y = x1
		C = np.matrix('1 0')
		# print C,type(C),x, type(x)
		y = C*x
		return y

	G.der = der
	G.out = out
	G.namestring = "Plant G"

	G.cstates = 2
	G.inargs = 1
	G.passargs = []

def K():

	def der(t,x,x_ref,xhat):
		return np.array([])
	
	def out(t,x,x_ref,xhat):
		# print x_ref[0,0], xhat
		assert len(x_ref) == 1
		# e = np.concatenate( (np.array(x_ref),np.array([0])), axis=0)-xhat
		e = np.matrix('1 ; 0')*x_ref-xhat
		Kmat = np.matrix('1 2')
		# u = Kmat*np.reshape(e,(2,1))
		u = Kmat*e
		# u = u[0,0]
		return u
		# return np.array([u[0,0]])

	K.der = der
	K.out = out
	K.namestring = "Controller K"

	K.inargs = 2
	K.cstates = 0
	K.passargs = [0,1]
	

def L():

	def der(t,x,u,y):
		return np.array([ x[1], -1.0*x[0] -2.0*x[1] + u[0] ])
	def out(t,x,u,y):
		L = np.matrix('0.5; 1.5')
		return 1.0*x+ 0.0*L*y
	
	L.der = der
	L.out = out
	L.namestring = "Observer L"

	L.inargs = 2
	L.cstates = 2
	L.passargs = [1]


def ES():
	#external source

	def der(t,x):
		return np.matrix([])
	def out(t,x):
		return np.matrix([np.sin(t)])
	
	ES.der = der
	ES.out = out
	ES.namestring = "external input"

	ES.inargs = 0
	ES.cstates = 0
	ES.passargs = []



def test1():
	# x_ref = np.matrix([[1.0], [2.0], [4.0]])
	

	# T = np.arange(0,1.0,0.01)
	# sys = (sll.integrator_multi(5,3),[])
	# x0 = ([0.0]*5*3,[])

	x_ref = np.matrix([[1.0]])

	A = np.matrix('0 1; -5 -0.01')
	B = np.matrix('0; 5')
	C = np.matrix('1 0')
	D = np.matrix('1')


	Trange = np.arange(0,10.0,0.001)
	sys = (sll.ss_LTI(A,B,C,D),[])
	x0 = ([0.0,0.0],[])

	M,x0 = sl.init_MDL(sys,x0,"this")

	x0 = np.transpose(np.matrix(x0))
	# Tad,Xad = ode.rungekutta4ad(M.der, x_ref, Trange, x0, e_tol=1e-5, min_dt=1e-4 )
	T,X = ode.rungekutta4(M.der, x_ref, Trange, x0)
	Y = [ M.out(t,x,x_ref) for t,x in zip(T,X) ]
	
	print T[-1], len(T)
	# print Tad[-1], len(Tad)
	print X[-1]
	print Y[-1]
	plt.figure(1)
	plt.subplot(211)
	plt.plot(T,[np.array(x)[:,0] for x in X])
	# plt.plot(Tad,[np.array(x)[:,0] for x in Xad])
	# plt.plot(T,[np.array(sl.go_deep(y[0]))[:,0] for y in Y] )

	plt.subplot(212)
	plt.plot(T)
	# plt.plot(Tad)

	plt.show()








def test2():
	fig, (ax1, ax2) = plt.subplots(2,1)
	states = ax1.plot( *([[]]*2) )
	outs = ax2.plot( *([[]]*2) )

	ax1.grid()

	def animate(T,X,*args):
		def hold():
			plt.show()

		

		animate.hold = hold
		ax1.figure.canvas.draw()
		ax2.figure.canvas.draw()

		this_x = np.array(X)
		out_probes = [0]
		if args:
			this_y = np.array([ [y.probe(p)[0,0] for p in out_probes] for y in args[0]])
			
		
		ax1.set_xlim(0, max(T))
		ax2.set_xlim(0, max(T))
		ax1.set_ylim(np.amin(this_x),np.amax(this_x))
		ax2.set_ylim(np.amin(this_y),np.amax(this_y))

		for i,l in enumerate(states):
			l.set_data(T, this_x[:,i] )
		for i,l in enumerate(outs):
			l.set_data(T,  this_y[:,i] )

		plt.pause(0.00001)


	T = np.arange(0,4.0,0.01)
	sys = (G,K,sll.gain(3.7),sll.add,sll.sgn,sll.fun_gen(A=3.0),sll.gain(5.0),sll.const(1.0),L,1,0)
	x0 = ([1.0, 0.6],[],[],[],[],[0.0,0.0],[],[],[])

	sys = (sll.sub,sll.int1,sll.fun_gen(A=0.0,omega=0.0,bias=1.0),sll.time)
	x0 = ([],[0.0],[],[],[])


	# sys = (G,K,x_ref,L,0,1)
	M,x0 = sl.init_MDL(sys,x0,"this")

	x0 = np.transpose(np.matrix(x0))
	T,X,Y = ode.rungekutta4ad(M.der, T, x0 , outcall=M.out, \
		adaptive=False, min_dt=1e-3, e_tol=1e-4, realtime=True, plotcall=animate)

	plt.show()


	
	
	print T[-1]
	print X[-1]
	# print Y[-1]

	# animate.hold()

def test4():


	T = np.arange(0,5.0,0.01)
	d_in = np.matrix([1.0])

	sys_1 = (G,K,sll.gain(3.7),[],L,1,0)
	x0_1 = ([1.0, 0.6],[],[],[],[0.0,0.0])

	M1,x0_1 = sl.init_MDL(sys_1, x0_1, "Model KGL1")

	sys_2 = (G,K,M1,[],L,1,0)
	x0_2 = ([1.5, 0.6],[],x0_1,[],[0.0,0.0],[],[])

	M2,x0_2 = sl.init_MDL(sys_2, x0_2, "Model KGL21")

	sys_12 = (M1,M2,D,ES,[])
	x0_12 =(x0_1,x0_2,[0.0,0.0],[],[],[])

	M,x0 = sl.init_MDL(sys_12 ,x0_12, "Model KGL1+KGL21-G")

	# M = sl.unpack_MDL(M)

	# M.out(0.0,x0,d_in)
	x0 = np.transpose(np.matrix(x0))
	T,X = ode.rungekutta4ad(M.der, d_in,  T, x0)
	Y = [ M.out(t,x,d_in) for t,x in zip(T,X) ]
	
	print T[-1]
	print X[-1]
	print Y[-1]
	print Y[-1].probe(1,2,3)
	plt.plot(T,[np.array(x)[:,0] for x in X] )
	plt.show()





if __name__ == '__main__':
	# print(sl.__doc__)
	test2()
