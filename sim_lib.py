#!/usr/bin/python

"""A simple python library for some common functions/blocks for sim_link
most functions are written and tested to support numpy matrices"""

import numpy as np



def int1():
	""" 1d, 1st order integrator"""
	def der(t,x,u):
		xdot = u
		return xdot
	def out(t,x,u):
		y = x[0:1]
		return y

	int1.der = der
	int1.out = out
	int1.namestring = "int1"

	int1.inargs = 1
	int1.cstates = 1
	int1.passargs = []


def add():
	""" add two signals out=a+b """
	# def der(t,x,a,b):
	# 	return x[0:0]
	def out(t,x,a,b):
		return a+b

	# add.der = der
	add.out = out
	add.namestring = "add"

	add.inargs = 2
	add.cstates = 0
	add.passargs = [0,1]


def sub():
	""" subtract two signals out=a+b """
	# def der(t,x,a,b):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(t,x,a,b):
		return a-b

	# sub.der = der
	sub.out = out
	sub.namestring = "sub"

	sub.inargs = 2
	sub.cstates = 0
	sub.passargs = [0,1]


def mul():
	""" multiply two signals out=a*b """
	# def der(t,x,a,b):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(t,x,a,b):
		return a*b

	# mul.der = der
	mul.out = out
	mul.namestring = "mul"

	mul.inargs = 2
	mul.cstates = 0
	mul.passargs = [0,1]


def div():
	""" divide two signals out=a/b """
	# def der(t,x,a,b):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(t,x,a,b):
		return a/b

	# div.der = der
	div.out = out
	div.namestring = "div"

	div.inargs = 2
	div.cstates = 0
	div.passargs = [0,1]


def inv():
	""" inverse a signal out=1/u """
	# def der(t,x,u):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(t,x,u):
		return 1.0/u

	# inv.der = der
	inv.out = out
	inv.namestring = "inv"

	inv.inargs = 1
	inv.cstates = 0
	inv.passargs = [0]



def sat():
	""" saturate by -1,1 """
	# def der(t,x,u):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(t,x,u):
		return u* (-1 < u)*(u < 1) + 1*(u >= 1) - 1*(u <= -1)

	# inv.der = der
	sat.out = out
	sat.namestring = "sat"

	sat.inargs = 1
	sat.cstates = 0
	sat.passargs = [0]

def sgn():
	""" signum function """
	# def der(t,x,u):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(t,x,u):
		return  1*(u > 0) - 1*(u < 0)

	# inv.der = der
	sgn.out = out
	sgn.namestring = "sgn"

	sgn.inargs = 1
	sgn.cstates = 0
	sgn.passargs = [0]

def time():
	""" output time as a signal """
	# def der(t,x,a,b):
	# 	return x[0:0]
	def out(t,x):
		return np.matrix([t])

	time.out = out
	time.namestring = "time"

	time.inargs = 0
	time.cstates = 0
	time.passargs = []


def sin():
	""" sin(u) """
	def out(t,x,u):
		return np.sin(u)

	sin.out = out
	sin.namestring = "sin"
	sin.inargs = 1
	sin.cstates = 0
	sin.passargs = [0]

def cos():
	""" cos(u) """
	def out(t,x,u):
		return np.cos(u)

	cos.out = out
	cos.namestring = "cos"
	cos.inargs = 1
	cos.cstates = 0
	cos.passargs = [0]

def exp():
	""" exp(u) """
	def out(t,x,u):
		return np.exp(u)

	exp.out = out
	exp.namestring = "exp"
	exp.inargs = 1
	exp.cstates = 0
	exp.passargs = [0]

def log():
	""" log(u) """
	def out(t,x,u):
		return np.log(u)

	log.out = out
	log.namestring = "log"
	log.inargs = 1
	log.cstates = 0
	log.passargs = [0]


class gain(object):
	""" gain block out = k*u """
	def __init__(self,k):
		self.k = k

		self.namestring = "gain (" + str(self.k) + ")"

		self.inargs = 1
		self.cstates = 0
		self.passargs = [0]

	# def der(self,t,x,u):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(self,t,x,u):
		y = self.k*u
		return y

class const(object):
	""" const block out = C """
	def __init__(self,C):
		self.C = C

		self.namestring = "const (" + str(self.C) + ")"

		self.inargs = 0
		self.cstates = 0
		self.passargs = []

	# def der(self,t,x):
	# 	xdot = x[0:0]
	# 	return xdot
	def out(self,t,x):
		return self.C


class power(object):
	"""docstring for power"""
	def __init__(self, n):
		# super(power, self).__init__()
		self.n = n

		self.namestring = "power (" + str(self.n) + ")"

		self.inargs = 1
		self.cstates = 0
		self.passargs = [0]

	def out(self,t,x,u):
		return np.power(u,self.n)


class selector(object):
	"""docstring for selector"""
	def __init__(self, switches):
		# super(selector, self).__init__()
		self.sw = switches

		self.namestring = "selector (" + str(self.sw) + ")"

		self.inargs = 1
		self.cstates = 0
		self.passargs = [0]

	def out(self,t,x,u):
		return u[self.sw]
		


class integrator_multi(object):
	"""integrator, multi dimensional, nth order"""
	def __init__(self, ord, u_dim):
		# super(integrator, self).__init__()

		self.ord = ord

		self.cstates = ord*u_dim
		self.inargs = 1
		self.passargs = [0] if self.ord == 0 else []

		self.namestring = "integrator_multi"

	def der(self,t,x,u):
		# print x
		dx = np.roll(x,-1,axis=0)
		assert len(u) == u_dim
		dx[self.ord-1::self.ord] = u
		return dx

	def out(self,t,x,u):
		# print x
		return x[0::self.ord]
		


class ss_LTI(object):
	"""
	LTI system described by A,B,C,D matrices
	normal usage 	ss_LTI(A,B,C)
					ss_LTI(A,B,C,D)

	"""
	def __init__(self,*args):
		
		self.A = args[0]
		self.B = args[1]
		self.C = args[2]
		if len(args)==4:
			self.D = args[3]
			nz = np.nonzero(self.D)
			if np.size(nz)==0:
				self.passargs = []
			else:
				self.passargs = [0]
		else:
			self.D = np.zeros((np.shape(self.C)[0],np.shape(self.B)[1]))
			self.passargs = []

		self.namestring = "ss_LTI"

		self.inargs = 1
		self.cstates = len(self.A)
		

	def der(self,t,x,u):
		xdot = self.A*x + self.B*u
		# print x, xdot
		return xdot
	def out(self,t,x,u):
		y = self.C*x + self.D*u
		return y

class fun_gen(object):
	"""y = A*cos(wt+phi)+ofs"""
	defaultdata = {"A": 1.0, "omega": 2.0*np.pi, "phi": 0.0, "bias": 0.0}
	def __init__(self,**kwargs):
		self.data = self.defaultdata.copy()
		self.data.update(kwargs)

		self.namestring = "fun_gen"

		self.inargs = 0
		self.cstates = 0
		self.passargs = []


	# def der(self,t,x):
	# 	return np.matrix([])
	def out(self,t,x):
		return np.matrix([self.data["A"]*np.cos(self.data["omega"]*t\
			+self.data["phi"]) + self.data["bias"]])

class joystick_input(object):
	""" joystick input for linux and windows """
	def __init__(self,*args):
		self.namestring = "joystick_input"

		self.inargs = 0
		self.cstates = 0
		self.passargs = []

		from inputs import get_gamepad
		self.get_gamepad = get_gamepad

		self.memory = 0.0

	# def der(self,t,x):
	# 	return np.matrix([])

	def out(self,t,x):
		events = self.get_gamepad()
		for event in events:
			# print(event.ev_type, event.code, event.state)
			if event.code == 'ABS_RZ':
				self.memory = -(1.0*event.state-128)/128
		joy_out = [self.memory]

		joy_out = np.matrix(joy_out)
		return joy_out


# def plot(T,X,*args,**kwargs):
# 	data = {"states": [], "probes": []}
# 	plot.statics = {"stateplots": [], "outplots": [], "plots_exist", False}
# 	data.update(kwargs)

# 	if plot

# 	for s in states:






	
	


