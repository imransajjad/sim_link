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

Step 1: Write functions representing models of dynamic systems and then create
MDL-like objects, for example:

---

def der(t,x,u,y):
    # x1dot = x2dot
    # x2dot = -1.0*x1 -2.0*x2 + u
    return np.array([ x[1], -1.0*x[0] -2.0*x[1] + u[0] ])
def out(t,x,u,y):
    # xhat = x + L*y
    return 0.5*np.array(x)+ np.array([0.5, 1.5])*np.array(y)

L = MDLBase(der, out, 2, [1], [0.0, 0.0], name="Observer L")
# 2 is number of inputs (u and y)
# [1] is indices of inputs that are required to compute output (y)
# [0.0, 0.0] is the initial state vector inside the system
# name="Observer L" is a human readable name to be used to identify the block

---


Step 2: Connect such systems together and evaluate

When providing a tuple to the system think of it as unpacking an expression:

G(K(diff(x_ref,L(K(...),G(...)))))  => (G,K,diff,[],L,1,0)

Note: There is no guarantee that an arbitrary expression like the one on the
right will work, but I've tried my best to make the ones on the left work for
all cases

elements in this tuple can only be:

function objects 	<=>	system blocks
empty lists []		<=>	inputs
ints e.g 0, 1		<=> references to other signals

The order of elements here is important!!!

For example, here is test code for the system:

          __     ------         ------
x_ref -->|+-|-->|  K   |------>|  G   |-------->
          --     ------   |     ------   |
          ^               V              |
          |             ------           |
           ------------|  L   |<---------
                        ------

If this is taken to be the usual plant, observer controller system,
the signals on the wires should be evaluated in the order G,L,K

---

import numpy as np
import ode_solvers as ode
import sim_link as sl


def test8():
    import matplotlib.pyplot as plt
    x_ref = [1.0]

    T = np.arange(0,10.0,0.01)
    sys = (G,K,diff,[],L,1,0)
    # sys = (G,K,diff,[],L,0,1) # try this, it'll cause an error

    M = sl.MDL(sys,x0,"this")
    x0 = np.transpose(np.matrix(M.x0))

    T,X = ode.rungekutta4ad(M.der, x_ref,T, x0 )
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
    MDL builds an execution table which der and out then use
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
        A(A(in)) => (A,A,[in]) => 1/s^2
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

import copy

verbose = False
def print_debug(x):
    if verbose:
        print(x)

class XList(object):
    """
    A list like data type that uses np arrays as its atomic objects and
    supports addition and scalar multiplication.
    Members are accessible like a regular list i.e. x[0], x[1] or if
    index_strs = ["id_0", "id_1"] is provided at construction,
    dictionary like access, x["id_0"], x["id_1"] also works
    """

    def __init__(self, init_list, index_strs=None):
        self._list = list(init_list)
        self._index_strs = index_strs
        if index_strs:
            self._index = { str_i:i for i,str_i in enumerate(index_strs)}

    def __add__(self, other):
        """
        only support adding scalar 0. otherwise only XList can be added to XList
        """
        if other == 0:
            return XList(self._list, self._index_strs)
        return XList([a + b if not (a is None or b is None) else None for a, b in zip(self._list, other)], index_strs=self._index_strs)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        """
        only support scalar multiplication
        """
        assert not isinstance(other, XList)
        return XList([a * other if not (a is None) else None for a in self._list], index_strs=self._index_strs)

    def __lmul__(self, other):
        return self.__mul__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        only support scalar division on right
        """
        assert not isinstance(other, XList)
        return XList([a / other if not (a is None) else None for a in self._list], index_strs=self._index_strs)

    def __getitem__(self, index):
        if type(index) == str:
            index = self._index[index]
        return self._list.__getitem__(index)

    def __setitem__(self, key, value):
        if type(key) == str:
            key = self._index[key]
        return self._list.__setitem__(key, value)

    def __repr__(self):
        return f"XList: keys: {self._index_strs.__repr__()} values: {self._list.__repr__()}"


class SignalType(object):
    """
    A small helper class for entries in MDL.signals_valid
    """

    UNKNOWN = 0  # (unknown)
    INPUT = 1  # is input to some function
    REF = 2  # is reference to something else
    OUT = 4  # output available

    @staticmethod
    def isvalid(value):
        return value & (SignalType.REF | SignalType.OUT)


class MDLBase(object):
    """
    MDLBase serves as the base for other model type objects
    """

    def __init__(self, der, out, n_inargs, passargs, x0, name="MDLBase"):
        self.name = name

        self.der = der  # the function to calculate the derivates of the state vector
        self.out = out  # the output of the system

        self.passargs = passargs  # the list of passargs is required
        self.inargs = n_inargs  # number of input arguments
        self.x0 = x0  # the initial state vector

    def get_x0(self):
        return self.x0

    def __repr__(self):
        return f"MDL:{self.name}"


class MDL(MDLBase):
    """
    MDL is an MDLBase type of object that combines other MDL type objects into a system
    """

    def __init__(self, sys, name=""):
        self.name = name
        self.der = self.mdl_der
        self.out = self.mdl_out

        print_debug(f"Initializing System: {self.name} with\n{sys}")

        # these values are used to only verify correctness
        self.signals_valid = [SignalType.UNKNOWN for i in sys]

        # these value are used in execution
        self.exec_table = []  # a table of the MDL type objects in sys only, ordered
        # each entry is (sys, inputs_map, states_map, output_map)
        self.exec_table_outrow = -1  # the location of the output in the exec table
        self.signal_reg = list(sys)  # a list of all the signals in the model

        self.passargs_map = []

        # keep building table till complete
        signals_valid = [SignalType.isvalid(i) for i in self.signals_valid]
        while not all(signals_valid):
            first_invalid_index = [i for (i, sig) in enumerate(signals_valid) if not sig][0]
            print_debug("\nbuilding from __init__")
            self.build_table(first_invalid_index)
            signals_valid = [SignalType.isvalid(i) for i in self.signals_valid]

        self.argmap = [i for i, r in enumerate(self.signal_reg) if isinstance(r, list)]

        self.inargs = len(self.argmap)
        self.passargs = [i for i, p in enumerate(self.argmap) if p in self.passargs_map]

        self.sys_list = [s.name if isinstance(s, MDLBase) else None for s in self.signal_reg]
        for i,j in enumerate(self.argmap):
            self.sys_list[j] = f"{self.name}_u_{i}"

        self.table_valid = self.verify_table()

    def eval_table(self, N):
        """
        Version of build table that does no writes and just tells the last entry
        that is required for computing Nth entry
        entry
        if sys is (G,K,diff,[],L,1,0,[]) and G.passargs = [1], then
        [eval_table(i) for i in range(0,8)] == [7,6,6,3,6,5,6,7]
        """
        obj = self.signal_reg[N]

        if isinstance(obj, MDLBase):
            j = N
            for i in range(0, obj.inargs):
                j += 1
                j = self.eval_table(j)
            return j
        elif isinstance(obj, int):
            return N
        elif isinstance(obj, list):
            return N

    def resolve_reference(self, N):
        obj = self.signal_reg[N]
        if isinstance(obj, int):
            # print(f"eval table returning {N}")
            return obj
        else:
            return N

    def build_table(self, N):
        """
        generate exec table recursively starting at the Nth signal in MDL.signal_reg
        return the index of MDL.signal_reg that was just evaluated or in case of
        a reference, the index of what was refered to
        """
        obj = self.signal_reg[N]

        argmap = []
        print_debug(f"build table {N}")

        if isinstance(obj, MDLBase) and not SignalType.isvalid(self.signals_valid[N]):
            # element is a function, need to add to exec_table
            print_debug(f"found function {obj} at {N}")
            self.signals_valid[N] |= SignalType.OUT  # optimistically say output will be available

            # then if obj has passargs, make sure where they come from is evaluated first

            j = N
            for i in range(0, obj.inargs):

                # i is arg number
                # j is next position where arg can be

                j += 1  # start searching from here
                argmap.append(self.resolve_reference(j))
                if i in obj.passargs:
                    self.build_table(j)
                    if self.exec_table_outrow < 0:
                        self.passargs_map.append(j)
                j = self.eval_table(j)

            new = (obj, argmap, N)
            if N == 0:
                print_debug("\tFound mdl output")
                self.exec_table_outrow = len(self.exec_table)
            self.exec_table.append(new)
            return j

        elif isinstance(obj, int):
            print_debug(f"found reference to {obj} at {N}")
            self.signals_valid[N] |= SignalType.REF
            return self.build_table(obj)
        elif isinstance(obj, list):
            self.signals_valid[N] |= SignalType.OUT
            print_debug(f"found constant {obj} at {N}")
            return N
        else:
            print_debug(f"found else {obj} at {N}")
            return N

    def verify_table(self):
        """
        Find algebraic loops in the system
        """
        p_inputs_table = [([argmap[i] for i in obj.passargs], obj) for obj, argmap, _ in self.exec_table]
        for row, (obj, _, output) in enumerate(self.exec_table):
            for i, (p_input, p_input_obj) in enumerate(p_inputs_table[0 : (row + 1)]):
                if output in p_input:
                    return f"Algebraic loop: output ({output}) of entry {row}:{obj} is required by entry {i}:{p_input_obj}"
        return "Valid"

    def get_x0(self):
        """
        Get the state of the system in a list of state vectors in the same order
        as the system tuple
        """
        state_list = [s.get_x0() if isinstance(s, MDLBase) else None for s in self.signal_reg]
        x0 = XList(state_list, index_strs=self.sys_list)
        return x0

    def mdl_out(self, t, x, *inputs, probes=False):
        """
        evaluate only till output is reached and return output
        """
        y = list(self.signal_reg)
        for i, arg_pos in enumerate(self.argmap):
            y[arg_pos] = inputs[i]
        for obj, inputs_map, output in self.exec_table:
            obj_u = [y[i] for i in inputs_map]
            y[output] = obj.out(t, x[output], *obj_u)
            if not probes and output == 0:
                return y[0]
        return XList(y, index_strs=self.sys_list)

    def mdl_der(self, t, x, *inputs):
        y = self.mdl_out(t, x, *inputs, probes=True)
        dx = copy.copy(x)
        for obj, inputs_map, output in self.exec_table:
            if not obj.der:
                continue
            obj_u = [y[i] for i in inputs_map]
            dx[output] = obj.der(t, x[output], *obj_u)
        return dx

    def table(self):
        string = ""
        string += (
            f"--------------\n"
            + f"System Execution Table ({self.name}) inargs: {self.inargs} "
            + f"passargs: {self.passargs}\n"
            # + f"state: {self.get_x0()}\n"
        )
        for obj, inputs, output in self.exec_table:
            i_list = [f"{inp}" + ("p" if i in obj.passargs else "") for i, inp in enumerate(inputs)]
            string += f"\t{obj.name}, ins:[{', '.join(i_list)}], out: {output}\n"
        string += f"Table valid: {self.table_valid}\n"
        string += "End System Execution Table (" + self.name + ")\n--------------"
        return string


def unpack_MDL(M):
    # if a model consists of submodels, this will 'unpack' it's execution table
    # until it only contains functions
    # does (should) not affect execution order or state variables

    def offset_row(sub_row, ofs):
        return [sub_row[0], [i + ofs for i in sub_row[1]], sub_row[-1] + ofs]

    def find_replace_inrow(sub_row, a, b):
        # replace a with b in row
        return [sub_row[0], [b if i == a else i for i in sub_row[1]], b if sub_row[-1] == a else sub_row[-1]]

    def reg_extend(M):
        Nreg = []
        Nreg += M.ETregister
        Nrows = list(M.ET)
        OFSes = []
        i = 0
        while i < len(Nreg):
            if isinstance(Nreg[i], MDL):

                # print "\n\n\n"
                # print i
                # print Nreg
                # print "\n"
                # for row in Nrows:
                # 	print row

                # find index of matching row in ET
                m = [j for j, row in enumerate(Nrows) if row[-1] == i][0]

                # print m

                sub_argmap = Nreg[i].argmap
                sub_reg = Nreg[i].ETregister
                # done extracting some data from subsystem

                ofs = len(Nreg) - 1

                OFSes.append([i, ofs, Nreg[i].inargs, Nreg[i].argmap])

                # start building replacement and additions to Nreg
                Nreg[i] = sub_reg[0]

                add_reg = [j + ofs if isinstance(j, int) else j for j in sub_reg[1:]]
                add_reg = [i if j == ofs else j for j in add_reg]

                arg_locs = [j for j, a in enumerate(add_reg) if isinstance(a, list)]

                assert len(arg_locs) == len(sub_argmap)

                # print "\n"
                # print Nreg

                # then replace rows in ET
                add_rows = [offset_row(row, ofs) for row in Nrows[m][0].ET]

                # replace output
                add_rows = [find_replace_inrow(row, ofs, i) for row in add_rows]
                # replace inputs
                for ins, sub_ins, reg_argi in zip(Nrows[m][1], sub_argmap, arg_locs):
                    add_reg[reg_argi] = ins
                    add_rows = [find_replace_inrow(row, sub_ins + ofs, ins) for row in add_rows]

                Nrows = Nrows[0:m] + add_rows + Nrows[m + 1 :]
                Nreg += add_reg
            i += 1

        return Nreg, Nrows, OFSes

    print("Unpacking System:", M.namestring)

    Nreg, Nrows, ofses = reg_extend(M)

    # for i in Nrows:
    # 	print i

    # print ofses
    # Nrows = list(M.ET)
    # outs_to_unpack = [i[0] for i in ofses]
    # for ofs in ofses:
    # 	m = [ i for i,j in enumerate(M.ET) if j[2]==ofs[0]][0]
    # 	# find index of matching row
    # 	add_rows = [offset_row(row,ofs[1]) for row in Nrows[m][0].ET]

    # 	# replace output
    # 	add_rows = [find_replace_inrow(row,ofs[1],ofs[0]) for row in add_rows]
    # 	# replace input
    # 	for ins,sub_ins in zip(Nrows[m][1],ofs[3]):
    # 		add_rows = [find_replace_inrow(row,sub_ins+ofs[1],ins) \
    # 				for row in add_rows]

    # 	Nrows = Nrows[0:m] + add_rows + Nrows[m+1:]

    # Mnew = MDL(M.namestring+"_unpacked")
    Mnew = copy.copy(M)
    Mnew.namestring += "_unpacked"

    Mnew.ET = Nrows
    Mnew.ETregister = Nreg
    Mnew.table_rep_stride()
    Mnew.verify_table()
    Mnew.ETvalid = [0] * len(Nreg)

    # Mnew.print_table()

    assert Mnew.ETflag[0]

    # print min_tuples
    # print Mnew.argmap

    # print "Done Unpacked System:", Mnew.namestring, "\n"

    return Mnew


def verify(A):
    """necessary test if A is a valid function block"""
    maybe_valid = True

    attr_names = ["out", "namestring", "inargs", "passargs", "cstates"]
    nec_attr = [hasattr(A, i) for i in attr_names]
    der_attr = hasattr(A, "der")

    if not all(nec_attr):
        maybe_valid = False
        print(A, "does not have:")
        # print nec_attr
        for name, a in zip(attr_names, nec_attr):
            if not a:
                print(name)

    if nec_attr[0] and der_attr:
        if not A.der.func_code.co_argcount == A.out.func_code.co_argcount:
            maybe_valid = False
            print(A, "out and der don't have equal number of arguments")
    if nec_attr[2] and nec_attr[3]:
        if not len(A.passargs) <= A.inargs:
            maybe_valid = False
            print(A, "passargs/inargs not right")
    return maybe_valid
