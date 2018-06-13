#!/usr/bin/python

"""
sim_link.py provides some functionality to connect and evaluate some
dynamical systems

Maintained on Github
https://github.com/imransajjad/sim_link

by
Imran Sajjad
imransajjad89@gmail.com


Basic workflow:

Step 1: Write functions representing models of dynamic systems, for example:

---

def L():

	def der(t,x,u,y):
		# x1dot = x2dot
		# x2dot = -1.0*x1 -2.0*x2 + u
		return np.array([ x[1], -1.0*x[0] -2.0*x[1] + u[0] ])
	def out(t,x,u,y):
		# xhat = x + L*y
		return 0.5*np.array(x)+ np.array([0.5, 1.5])*np.array(y)
	
	L.der = der
	L.out = out

	L.inargs = 2 # number of continous states
	L.cstates = 2 # number of input arguments in addition to t and x
	L.passargs = [1] # indices of input arguments actually used in L.out
					# in this case, y is a passthrough argument

---

L() should be callable (to initialize) and L.der and L.out should have the same
input arguments.



Step 2: Connect such systems together and evaluate

When providing a tuple to the system think of it as unpacking an expression:

G(K(x_ref,L(K(...),G(...))))  => (G,K,x_ref,L,1,0)

I cannot emphasize this enough. There is no guarantee that an arbitrary 
expression like the one on the right will work, but I've tried my best
to make the ones on the left work for all cases

elements in this tuple can only be:

function objects 	<=>	system blocks
empty lists []		<=>	inputs
ints e.g 0, 1		<=> references to other signals

The order of elements here is important!!!

For example, here is test code for the system:

           ------         ------
x_ref --->|  K   |------>|  G   |-------->
          ------    |     ------   |
             ^      V              |
             |    ------           |
              ---|  L   |<---------
                  ------

If this is taken to be the usual plant, observer controller system,
the signals on the wires should be evaluated in the order G,L,K

---

import numpy as np
import ode_solvers as ode
import sim_link as sl
import matplotlib.pyplot as plt


def test2():
	x_ref = [1.0]
	
	T = np.arange(0,10.0,0.01)
	sys = (G,K,[],L,1,0)
	x0 = ([1.0, 0.6],[],[],[0.0,0.0],[],[],[])
	# sys = (G,K,[],L,0,1) # try this, it'll cause an error

	M = sl.MDL(sys,x0,"this")
	x0 = np.transpose(np.matrix(M.x0))

	T,X = ode.rungekutta4(M.der, x_ref,T, x0 )
	Y = [ M.out(t,x,x_ref) for t,x in zip(T,X) ]
	
	print T[-1]
	print X[-1]
	print Y[-1]
	plt.plot(T,[np.array(x)[:,0] for x in X])
	plt.show()

---

Basic ideas:
1. t,x are global (in a sense)
    so construct a der(t,x,*args) and out(t,x,*args) by connecting individual
    functions, and 'extract' hidden states from subsystems into the global 
    state vector

2. evaluation order has to be correct for out
    no algebraic loops are allowed,
    so MDL builds an execution table which der and out then use
    order of der execution doesn't really matter as long as out has been 
    evaluated
    passargs are necessary to eliminate false positives in algebraic 
    loops

3. This is kind of like a Moore State machine
    The results of all systems are stored in a register which is then used
    by subsequent systems, in the proper order

4. Need so solve the differential equation formed by der and out:
    xdot = der(t,x,*args)
    y = out(t,x,*args)
    You can use someone else's solver, write your own, or use the ode_solvers 
    module

5. The MDL function is also similar to G, L etc
    It has MDL(), MDL.out(t,x), MDL.der(t,x) methods, which are used in the same
    way
    In addition, it has the static variables MDL.ET and MDL.ETflag[0] that contain
    references to the original functions, execution order and valid flags, that
    together form the execution table.

6. Multiple instances of the same function (block) are allowed:
    Since the state is not stored inside the functions themselves, 'two or more'
    of the same system will act as two separate systems and their states and 
    results will be stored separately, for example if A(s) = 1/s
        A(A(in)) <=> (A,A,[in]) => 1/s^2
    If you want to use the result of a function twice, use references (ints)

7. Since MDL objects also act as dynamical system functions
	They can also be connected together
	sys_1 = (G,K,[],L,1,0)
	M1,x0_1 = sl.init_MDL(sys_1, x0_1, "Model KGL1")
	
	sys_2 = (G,K,[],L,1,0)
	M2,x0_2 = sl.init_MDL(sys_2, x0_2, "Model KGL2")

	sys_12 = (D,M1,M2,D,ES,[],[])
	M,x0 = sl.init_MDL(sys_12 ,x0_12, "Model KGL1+KGL2-G")

	is a valid operation

8. MDL objects can be unpacked
	If a MDL object has references to execute other MDL objects, the 
	unpack_MDL() method will replace the row containing an MDL object
	with that MDL object's rows. This might be useful, but is not recommended
	because while the execution order and state vector order are maintained,
	the output vector order is not guaranteed to be the same.
	I have assumed that the original order cannot be extracted from a 'packed'
	system, but I may be wrong.

9. External Inputs are also allowed
	Simply define a function with no states or inputs, but as an output, you
	can make any type of system calls you want. The ode solvers will have to
	be run in realtime for this.



Future ideas:
realtime option in solver
fix adaptive solver
add discrete states and mixed solver
I don't want to use isinstance and hasattr so much

"""

import copy, traceback


class out_list(list):
	"""docstring for out_list"""
	def getitemnormally(self,index):
		return super(out_list,self).__getitem__(index)
	def getfirstitem(self, index):
		R = self.getitemnormally(index)
		if isinstance(R,out_list):
			R = R.getfirstitem(0)
		return R

	def __getitem__(self, index):
		return self.getfirstitem(index)
	def probe(self, *indices):
		if len(indices) > 1:
			return self.getitemnormally(indices[0]).probe(*indices[1:])
		else:
			return self.getfirstitem(indices[0])
		# the altered getitem method is experimental, if in doubt use
		# the go_deep() function
	def probe_s(self, *probe_indices):
		return [self.probe(*i) for i in probe_indices]
	

class MDL(object):
	"""docstring for MDL"""
	def __init__(self,sys,x0_in,namestring):
		# super(MDL, self).__init__()
		self.namestring = namestring


		print "Initializing System:", namestring
	
		for i in sys:
			if callable(i):
				try:
					i()
				except Exception as e:
					print "Could not initialize System:", i
					# traceback.print_exc()
					pass

		self.ETregister = out_list(sys)

		# evaluate test and remap initial conditions
		self.init_table() # build table

		assert self.ETflag[0]

		self.x0 = []
		for row in self.ET:
			self.x0 += x0_in[row[2]]

		
		self.print_table()
		# print "x0: ", self.x0
		# print "Done Initializing System:", namestring, "\n\n\n"

	
	def init_table(self):

		
		self.ETflag = [False]
		self.ET = []
		self.ETvalid = [0]*len(self.ETregister)

		# keep building table till complete
		while not all( [i&6 for i in self.ETvalid] ):
			# print "build table pass"
			self.build_table([i for (i,T) in enumerate(self.ETvalid) \
						if not T&6][0],True)
		
		self.ETflag[0] = self.verify_table()
		self.argmap = [i for i,r in enumerate(self.ETregister) \
						if isinstance(r,list)]
		
		self.inargs = len(self.argmap)
		self.cstates = sum([row[0].cstates for row in self.ET ])
		self.passargs = [i for i, p in enumerate([self.isPassArg(i,False) \
						for i in self.argmap]) if p]
		
		return self.ETflag[0]


	def build_table(self,N,write):
		# generate exec table recursively
		args = self.ETregister

		# print "\n", self.ETvalid
		exec_is = []

		# entries in valid table are
		# 0 (unevaluated)
		# 1 is input to some function
		# 2 is reference to something else
		# 4 output available

		verbose = False
		if verbose: print self.ETvalid

		
		if write and hasattr(args[N], 'out') and not self.ETvalid[N]&4:
			# base element is a function, need to add to table
			if verbose: print N, args[N], "found function"
			
			self.ETvalid[N] |= 4 # optimistically say output will be available
			for i in range(0,args[N].inargs):
				
				# print i in args[N].passargs

				# i is arg number
				# j is next position where arg can be

				j = N+1+i
				while self.ETvalid[j]&1:
					j += 1
					if j == len(args):
						print "not enough args:", N, "needs arg", j
						print "but len(args) is", len(args)
				self.ETvalid[j] |= 1
				j = self.build_table(j,i in args[N].passargs)
				
				exec_is.append(j)

			# print exec_is
			
			# print "writing: ", args[N], exec_is, N, "\n"
			self.ET.append( [args[N], exec_is, N ]  )
			
			return N

		elif isinstance(args[N], int):
			if verbose: print N, args[N], "found reference to", args[N]
			self.ETvalid[N] |= 2

			return self.build_table(args[N],write)
		elif isinstance(args[N], list):
			self.ETvalid[N] |= 4
			if verbose:  print N, args[N], "found constant", args[N]
			return N
		else:
			if verbose:  print N, args[N], "found else", args[N]
			return N


	def verify_table(self):

		try:
			
			inputs = [[row[1][i] for i in row[0].passargs] for row in self.ET ]
			# find row number containing
			for k,cur_row in enumerate(self.ET):
				J = [ i for i,ins in enumerate(inputs) if cur_row[2] in ins]
				# print "J: ",J
				if J:
					j = J[0]
					if j <= k:
						print "row", j, "invalid because of row", k , "...\n"
						raise ValueError('System Data invalid')
						
						
						
		except Exception as e:
			self.print_table()
			print "verification failed"
			traceback.print_exc()

		return True

	def isPassArg(self,argi,deep):
		# first find row(s)
		if deep:
			next_args = [ row[2] for row in self.ET \
					if argi in [ row[1][i] for i in row[0].passargs ]  ]
			
			# if no next arg needs argi
			if not next_args:
				return False
			elif any([n==0 for n in next_args]):
				return True
			else:
				return any([self.isPassArg(n) for n in next_args])
		else:
			return [ argi for row in self.ET if argi in [ row[1][i] \
					for i in row[0].passargs ]  ]

	def out(self,t,x,*inputs):
		R = out_list(self.ETregister) # y

		if not len(inputs) == self.inargs:
			print "I ", self, " need ", self.inargs, " but received:"
			print inputs
			raise AssertionError()
		for i in range(0,self.inargs):
			R[self.argmap[i]] = inputs[i]


		x_stride = 0
		for row in self.ET:
			try:
				# all R's should be first output signals, dive inside till true
				# ins = go_deep([R[i] for i in row[1]])
				ins = [R.getfirstitem(i) for i in row[1]]
				
				R[row[2]] = row[0].out(t,x[x_stride:x_stride+row[0].cstates], *ins )
				x_stride += row[0].cstates
			except Exception as e:
				traceback.print_exc()
				print "Error in out from this model, row:"
				print "I,", self.namestring, ", called"
				print row
				print "with states", x[x_stride:x_stride+row[0].cstates]
				print "with inputs", ins
				raise e

		return R

	def der(self,t,x,*inputs):
		R = self.out(t,x,*inputs)
		
		# # all R's should be first output signals, dive inside till true
		# for i,r in enumerate(R):
		# 	if not isinstance(R[i],int) and not hasattr(R[i],'out'):
		# 		while not isinstance(R[i][0],float):
		# 			R[i] = R[i][0]

		Dx = copy.copy(x)

		x_stride = 0
		for row in self.ET:
			try:
				# print row
				# print t,x, x_stride,x_stride+row[0].cstates, inputs

				# all R's should be first output signals, dive inside till true
				# ins = go_deep([R[i] for i in row[1]])
				if row[0].cstates:
					ins = [R.getfirstitem(i) for i in row[1]]

					dx = row[0].der(t,x[x_stride:x_stride+row[0].cstates], *ins )
					for xi in range(x_stride,x_stride+row[0].cstates):
						Dx[xi] = dx[xi-x_stride]
					x_stride += row[0].cstates
			except Exception as e:
				traceback.print_exc()
				print "Error in der from this model, row:"
				print self.namestring, row
				print [R[i] for i in row[1]]
				raise e

		return Dx

	def print_probes(self):
		def pp(self,ind):
			for i,r in enumerate(self.ETregister):
				print "%s%d %s" % (ind, i, ( r.namestring if hasattr(r,'out')  else "" ))
				if isinstance(r,MDL):
					pp(r,"   " + ind)
		print "Valid probe values:"
		pp(self,"")



	def print_table(self):
		print "--------------\nSystem Execution Table ("+self.namestring+ ")"
		for row in self.ET:
			i_string = [ str(i_s) for i,i_s in enumerate(row[1])]
			if not isinstance(row[0],list):
				for i in row[0].passargs:
					i_string[i] += "p"
			print "	", row[0].namestring, ", ins:" , i_string, ", outs:", row[2]
		print "Table valid: ", self.ETflag[0]
		self.print_register()
		if self.ETflag[0]:
			self.print_addinfo()
			self.print_probes()
		
		print "\nEnd System Execution Table ("+ self.namestring +")\n--------------\n"

	def print_register(self):
		print "\nSystem signal register:"
		print "	", [R.namestring if hasattr(R,'namestring') \
					else R for R in self.ETregister]

	def print_addinfo(self):
		print "\nSystem I/O:"
		i_string = [str(a)+"p" if i in self.passargs else str(a) \
						for i,a in enumerate(self.argmap)]
		print "Inputs", i_string
		print "Number of cstates", self.cstates, "\n"


def go_deep(ins):
	for i,r in enumerate(ins):
		if not isinstance(ins[i],int) and not hasattr(ins[i],'out'):
			while isinstance(ins[i],out_list):
				# print ins[i], ins[i][0]
				ins[i] = ins[i][0]
	return ins

		
# def init_MDL(sys,x0_in, namestring):

# 	print "Initializing System:", namestring, "\n"
	
# 	for i in sys:
# 		if callable(i):
# 			try:
# 				i()
# 			except Exception as e:
# 				print "Could not initialize System:", i
# 				# traceback.print_exc()
# 				pass

# 	# evaluate test and remap initial conditions
# 	M = MDL(namestring) # init and
# 	M.init_table(sys) # build table

# 	assert M.ETflag[0]

# 	x0 = []
# 	for row in M.ET:
# 		x0 += x0_in[row[2]]

	
# 	M.print_table()
# 	print "x0: ", x0
# 	print "Done Initializing System:", namestring, "\n\n\n"

	return M,x0

def unpack_MDL(M):
	# if a model consists of submodels, this will 'unpack' it's execution table
	# until it only contains functions
	# does (should) not affect execution order or state variables

	def offset_row(sub_row,ofs):
		return [sub_row[0],[i+ofs for i in sub_row[1]], sub_row[2]+ofs]

	def find_replace_inrow(sub_row,a,b):
		# replace a with b in row
		return [sub_row[0],[b if i==a else i for i in sub_row[1]], b \
						if sub_row[2]==a else sub_row[2]]


	def reg_extend(M):
		Nreg = []
		Nreg += M.ETregister
		Nrows = list(M.ET)
		OFSes = []
		i = 0
		while i < len(Nreg):
			if isinstance(Nreg[i],MDL):

				# print "\n\n\n"
				# print i
				# print Nreg
				# print "\n"
				# for row in Nrows:
				# 	print row

				# find index of matching row in ET
				m = [ j for j,row in enumerate(Nrows) if row[2]==i][0] 

				# print m
				
				sub_argmap = Nreg[i].argmap
				sub_reg = Nreg[i].ETregister
				# done extracting some data from subsystem

				ofs = len(Nreg)-1

				OFSes.append( [i,ofs,Nreg[i].inargs,Nreg[i].argmap] )

				# start building replacement and additions to Nreg
				Nreg[i] = sub_reg[0]
				
				add_reg = [j+ofs if isinstance(j,int) else j for j in sub_reg[1:]]
				add_reg = [i if j==ofs else j for j in add_reg]
				
				arg_locs = [j for j,a in enumerate(add_reg) if isinstance(a,list)]

				assert len( arg_locs ) == len( sub_argmap )



				# print "\n"
				# print Nreg

				# then replace rows in ET
				add_rows = [offset_row(row,ofs) for row in Nrows[m][0].ET]

				# replace output
				add_rows = [find_replace_inrow(row,ofs,i) for row in add_rows]
				# replace inputs
				for ins, sub_ins, reg_argi in zip( Nrows[m][1], sub_argmap, arg_locs ):
					add_reg[reg_argi] = ins
					add_rows = [find_replace_inrow(row,sub_ins+ofs,ins) for row in add_rows]

				Nrows = Nrows[0:m] + add_rows + Nrows[m+1:]
				Nreg += add_reg
			i += 1

		return Nreg, Nrows, OFSes

	print "Unpacking System:", M.namestring

	Nreg, Nrows, ofses = reg_extend(M)

	# for i in Nrows:
	# 	print i

	# print ofses
	# Nrows = list(M.ET)
	# outs_to_unpack = [i[0] for i in ofses]
	# for ofs in ofses:
	# 	m = [ i for i,j in enumerate(M.ET) if j[2]==ofs[0]][0]
	#	# find index of matching row
	# 	add_rows = [offset_row(row,ofs[1]) for row in Nrows[m][0].ET]

	# 	# replace output
	# 	add_rows = [find_replace_inrow(row,ofs[1],ofs[0]) for row in add_rows]
	# 	# replace input
	# 	for ins,sub_ins in zip(Nrows[m][1],ofs[3]):
	# 		add_rows = [find_replace_inrow(row,sub_ins+ofs[1],ins) \
	#				for row in add_rows]

	# 	Nrows = Nrows[0:m] + add_rows + Nrows[m+1:]

	# Mnew = MDL(M.namestring+"_unpacked")
	Mnew = copy.copy(M)
	Mnew.namestring += "_unpacked"

	Mnew.ET = Nrows
	Mnew.ETregister = Nreg
	Mnew.ETflag[0] = Mnew.verify_table()
	Mnew.ETvalid = [0]*len(Nreg)

	Mnew.print_table()

	assert Mnew.ETflag[0]

	# print min_tuples
	# print Mnew.argmap

	# print "Done Unpacked System:", Mnew.namestring, "\n"

	return Mnew

def verify(A):
	"""necessary test if A is a valid function block"""
	maybe_valid = True

	attr_names = ['out','namestring','inargs','passargs','cstates']
	nec_attr = [hasattr(A,i) for i in attr_names]
	der_attr = hasattr(A,'der')

	if not all(nec_attr):
		maybe_valid = False
		print A, "does not have:"
		# print nec_attr
		for name,a in zip(attr_names,nec_attr):
			if not a: print name

	if nec_attr[0] and der_attr:
		if not A.der.func_code.co_argcount == A.out.func_code.co_argcount:
			maybe_valid = False
			print A, "out and der don't have equal number of arguments"
	if nec_attr[2] and nec_attr[3]:
		if not len(A.passargs) <= A.inargs:
			maybe_valid = False
			print A, "passargs/inargs not right"
	return maybe_valid


