import numpy as np
import ode_solvers as ode
import sim_link as sl
import sim_lib as sll
import time

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
	import matplotlib.pyplot as plt

	T = np.arange(0,20.0,0.01)
	sys = (G,K,sll.gain(3.7),sll.add,sll.fun_gen(A=3.0,omega=2),sll.add,sll.gain(1.0),sll.function(np.random.randn,1),sll.step(5.0),L,1,0)
	x0 = ([1.0, 0.6],[],[],[],[],[],[],[0.0,0.0],[],[])

	# sys = (sll.sub,sll.int1,sll.fun_gen(A=1.0,omega=2*np.pi*2.0,bias=1.0),sll.time)
	# x0 = ([],[0.0],[],[],[])


	# sys = (G,K,x_ref,L,0,1)
	M = sl.MDL(sys,x0,"this")
	M.print_table()


	PW = sll.plot_window([0,1], [9,0], [9,1], plot_separate=False )
	# PW2 = sll.plot_window([0], [5], [2], [0], plot_separate=False, active_draw=True)
	fig, ax = PW.return_axes()
	ax[1].set_xlabel("no time?")

	x0 = np.transpose(np.matrix(M.x0))
	T,X,Y = ode.rungekutta4ad(M.der, T, x0 , outcall=M.all_out, \
		adaptive=False, min_dt=5e-3, e_tol=1e-4, realtime=True, plotcalls=[PW.animate])

	PW.show()
	plt.plot(T)
	plt.show()
	
	print T[-1]
	print X[-1]
	# print Y[-1]


def test4():
	import matplotlib.pyplot as plt


	T = np.arange(0,5.0,0.01)
	d_in = np.matrix([1.0])

	sys_1 = (G,K,sll.gain(3.7),[],L,1,0)
	x0_1 = ([1.0, 0.6],[],[],[],[0.0,0.0])

	M1 = sl.MDL(sys_1, x0_1, "Model KGL1")
	# M1.print_table()

	sys_2 = (G,K,M1,[],L,1,0)
	x0_2 = ([1.5, 0.6],[],M1.x0,[],[0.0,0.0],[],[])

	M2 = sl.MDL(sys_2, x0_2, "Model KGL21")
	# M2.print_table()

	sys_12 = (M1,M2,D,ES,[])
	x0_12 =(M1.x0,M2.x0,[0.0,0.0],[],[],[])

	M = sl.MDL(sys_12 ,x0_12, "Model KGL1+KGL21-G")
	assert sl.verify(M)
	M.print_table()

	# M = sl.unpack_MDL(M)
	# assert sl.verify(M)
	# M.print_table()

	# M.out(0.0,x0,d_in)
	x0 = np.transpose(np.matrix(M.x0))
	T,X = ode.rungekutta4ad(M.der, d_in,  T, x0)
	Y = [ M.out(t,x,d_in) for t,x in zip(T,X) ]
	
	print T[-1]
	print X[-1]
	print Y[-1]
	print Y[-1].probe_s([0],[2],[3])
	plt.plot(T,[np.array(x)[:,0] for x in X] )
	plt.show()


def invalid_function():
	def der(t,x):
		return 0
	invalid_function.der = der

	invalid_function.cstates = 0

def test0():
	sys = (sll.gain(3.4),[])
	x0 = ([],)


	G1 = sl.MDL(sys,x0,'sys_gain1')


	M = sl.MDL( (G1, sll.const(1.5)),([],[]),'sys1')
	M = sl.unpack_MDL(M)
	sys[0].k = 1.0


	# x0 = np.transpose(np.matrix(M.x0))
	t = np.arange(0,1,0.01)

	T,X = ode.rungekutta4ad(M.der,t,x0)
	Y = [M.out(t,x) for t,x in zip(T,X)]
	print T[-1]
	print X[-1]
	print Y[-1]



def test8():
	import matplotlib.pyplot as plt
	x_ref = [1.0]
	
	T = np.arange(0,10.0,0.01)
	sys = (G,K,[],L,1,0)
	x0 = ([1.0, 0.6],[],[],[0.0,0.0],[],[],[])
	# sys = (G,K,[],L,0,1) # try this, it'll cause an error

	M = sl.MDL(sys,x0,"this")
	x0 = np.transpose(np.matrix(M.x0))

	T,X = ode.rungekutta4ad(M.der, x_ref,T, x0 )
	Y = [ M.out(t,x,x_ref).probe_s([0],[1]) for t,x in zip(T,X) ]
	
	print T[-1]
	print X[-1]
	print Y[-1]
	plt.plot(T,[np.array(x)[:,0] for x in X])
	plt.show()


def test7():

	T = np.arange(-0.0,10.0,0.01)

	sys1 = (G,K,[],L,1,0,sll.int1,2)
	x01 = ([1.0, 0.6],[],[],[0.0,0.0],[],[],[0.0])
	M1 = sl.MDL(sys1,x01,'sys1')
	M1.print_table()

	sys2 = (G,K,[],L,1,0)
	x02 = ([0.0, 0.0],[],[],[0.0,0.0],[],[],[])
	M2 = sl.MDL(sys2,x02,'sys2')
	M2.print_table()

	x_in = np.matrix('8.0')

	sys = (sll.gain(5.2),M1,sll.sub,[],M2,0)
	x0 = ([],M1.x0,[],[],M2.x0,[])

	M = sl.MDL(sys, x0, 'false algebraic loop')
	M.print_table()
	# M = sl.unpack_MDL(M)
	# M.print_table()


	PW = sll.plot_window([0,2,4],[0], [1,0], plot_separate=True)
	fig, axes = PW.return_axes()
	axes[1].legend(['out','sys1-G'])
	
	x0 = np.transpose(np.matrix(M.x0))
	# y0 = M.all_out(T[0], x0, x_in)
	# print y0

	# print [ y0[i] for i in range(0,len(y0))]
	# dx0 = M.der(T[0], x0, x_in)
	# print dx0

	
	T,X,Y = ode.rungekutta4ad(M.der, x_in, T, x0 , outcall=M.all_out, \
		adaptive=False, realtime=False, plotcalls=[PW.animate], plottime=0.01)

	PW.show()
	
	print T[-1]
	print X[-1]
	print Y[-1].probe_s([0],[1,0],[1,1])








if __name__ == '__main__':
	# print(sl.__doc__)
	test2()
