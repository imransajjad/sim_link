#!/usr/bin/python

import numpy as np


# all signals are taken to be numpy.matrix objects

def int1():

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

	def der(t,x,a,b):
		return x[0:0]
	def out(t,x,a,b):
		y = a+b
		return y

	add.der = der
	add.out = out
	add.namestring = "add"

	add.inargs = 2
	add.cstates = 0
	add.passargs = [0,1]


def sub():

	def der(t,x,a,b):
		xdot = x[0:0]
		return xdot
	def out(t,x,a,b):
		y = a-b
		return y

	sub.der = der
	sub.out = out
	sub.namestring = "sub"

	sub.inargs = 2
	sub.cstates = 0
	sub.passargs = [0,1]


class gain(object):
	def __init__(self,k):
		self.k = k

		# self.der = der
		# self.out = out
		self.namestring = "gain (" + str(self.k) + ")"

		self.inargs = 1
		self.cstates = 0
		self.passargs = [0]

	def der(self,t,x,u):
		xdot = x[0:0]
		return xdot
	def out(self,t,x,u):
		y = self.k*u
		return y


class integrator_multi(object):
	"""docstring for integrator"""
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
		dx[self.ord-1::self.ord] = u
		return dx

	def out(self,t,x,u):
		# print x
		return x[0::self.ord]
		


class ss_LTI(object):
	def __init__(self,*args):
		# normal usage 	ss_LTI(A,B,C)
		#				ss_LTI(A,B,C,D)
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



