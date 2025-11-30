#!/usr/bin/python

"""A simple python library for some common functions/blocks for sim_link
most functions are written and tested to support numpy matrices"""

import numpy as np
import scipy

import sim_link.sim_link as sl


# @diagram_block(label="Gain", style="triangle;whiteSpace=wrap;html=1;", geometry=(-60, 150, 60, 80))
# class Gain(sl.MDLBase):
#     """gain block out = k*u"""

#     def __init__(self, k, name="gain"):
#         super().__init__(None, self.out, 1, [0], None, name=name)
#         self.k = k

#     def out(self, t, x, u):
#         y = self.k * u
#         return y

def integrator(name="integrator", x0=None):
    def der(t, x, u):
        xdot = u
        return xdot

    def out(t, x, u):
        y = x
        return y

    return sl.MDLBase(der, out, 1, [], x0, name=name)


def add(name="add"):
    """add two signals out=a+b"""

    def out(t, x, a, b):
        return a + b

    return sl.MDLBase(None, out, 2, [0, 1], None, name=name)


def add3(name="add"):
    """add three signals out=a+b+c"""

    def out(t, x, a, b, c):
        return a + b + c

    return sl.MDLBase(None, out, 3, [0, 1, 2], None, name=name)


def sub(name="sub"):
    """subtract two signals out=a-b"""

    def out(t, x, a, b):
        return a - b

    return sl.MDLBase(None, out, 2, [0, 1], None, name=name)


def mul(name="mul"):
    """multiply two signals out=a*b"""

    def out(t, x, a, b):
        return a * b

    return sl.MDLBase(None, out, 2, [0, 1], None, name=name)


def div(name="div"):
    """divide two signals out=a/b"""

    def out(t, x, a, b):
        return a / b

    return sl.MDLBase(None, out, 2, [0, 1], None, name=name)


def inv(name="inv"):
    """inverse a signal out=1/u"""

    def out(t, x, u):
        return 1.0 / u

    return sl.MDLBase(None, out, 1, [0], None, name=name)


def sat(name="sat"):
    """saturate by -1,1"""

    def out(t, x, u):
        return u * (-1 < u) * (u < 1) + 1 * (u >= 1) - 1 * (u <= -1)

    return sl.MDLBase(None, out, 1, [0], None, name=name)


def sgn(name="sgn"):
    """signum function"""

    def out(t, x, u):
        return 1 * (u > 0) - 1 * (u < 0)

    return sl.MDLBase(None, out, 1, [0], None, name=name)


def time(name="time"):
    """multiply input by time"""

    def out(t, x, u):
        return u * t

    return sl.MDLBase(None, out, 1, [0], None, name=name)


def sine(name="sine"):
    """take sine of input"""

    def out(t, x, u):
        return np.sin(u)

    return sl.MDLBase(None, out, 1, [0], None, name=name)


def cosine(name="cosine"):
    """take cosine of input"""

    def out(t, x, u):
        return np.cos(u)

    return sl.MDLBase(None, out, 1, [0], None, name=name)


def exp(name="exp"):
    """take exp of input"""

    def out(t, x, u):
        return np.exp(u)

    return sl.MDLBase(None, out, 1, [0], None, name=name)


def log(name="log"):
    """take log of input"""

    def out(t, x, u):
        return np.log(u)

    return sl.MDLBase(None, out, 1, [0], None, name=name)


class Gain(sl.MDLBase):
    """gain block out = k*u"""

    def __init__(self, k, name="gain"):
        super().__init__(None, self.out, 1, [0], None, name=name)
        self.k = k

    def out(self, t, x, u):
        y = self.k * u
        return y


class Observer(sl.MDLBase):
    """observer for an LTI State Space System"""

    def __init__(self, A, B, C, L, name="Obs"):
        super().__init__(self.der, self.out, 1, [0], None, name=name)
        self.A = A
        self.B = B
        self.C = C
        self.L = L

    def der(self, t, x, u, y):
        return self.A * x + self.B * u

    def out(self, t, x, u, y):
        return x + self.L * (y - self.C * x)


class StateSpaceLTI(sl.MDLBase):
    """an SS.LTI function defined by A,B,C and D matrices"""

    def __init__(self, A, B, C, D=None, x0=None, name="SS_LTI"):

        if x0 is None:
            x0 = np.array([[0]] * A.shape[1])
            print("x0 ", x0)
        super().__init__(self.der, self.out, 1, [0] if D else [], x0, name=name)
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def der(self, t, x, u):
        return np.matmul(self.A,x) + self.B*u

    def out(self, t, x, u):
        if self.D:
            return np.matmul(self.C,x) + np.matmul(self.D,u)
        else:
            return np.matmul(self.C,x)


class TransferFunction(StateSpaceLTI):
    """an SS.LTI function defined by a transfer function"""

    def __init__(self, num, den, form="control_canonical", name="TF"):
        A, B, C, D = scipy.signal.tf2ss(num, den)
        super().__init__(A, B, C, D=(D if len(num) == len(den) else None), name=name)


class Step(sl.MDLBase):
    """step block out = (t>ts)*in"""

    def __init__(self, ts, const=1.0, name="step"):
        super().__init__(None, self.out, 0, [], None, name=name)
        self.ts = ts
        self.const = const

    def out(self, t, x):
        return (t > self.ts) * self.const


class Const(sl.MDLBase):
    """Const output"""

    def __init__(self, const=1.0, name="const"):
        super().__init__(None, self.out, 0, [], None, name=name)
        self.const = const

    def out(self, t, x):
        return self.const


class Pow(sl.MDLBase):
    """Raise to a constant power"""

    def __init__(self, power, name="const"):
        super().__init__(None, self.out, 1, [0], None, name=name)
        self.power = power

    def out(self, t, x, u):
        return u**self.power


class SignalGen(sl.MDLBase):
    """y = A*cos(wt+phi)+ofs"""

    def __init__(self, amplitude, omega, phase=0, offset=0, name="signal_gen"):
        super().__init__(None, self.out, 0, [], None, name=name)
        self.A = amplitude
        self.omega = omega
        self.phase = phase
        self.offset = offset

    def out(self, t, x):
        return self.A * np.cos(self.omega * t + self.phase) + self.offset


class RandomGen(sl.MDLBase):
    """y = A*randn(args)"""

    def __init__(self, amplitude, randfn, randargs, name="rand_gen"):
        super().__init__(None, self.out, 0, [], None, name=name)
        self.A = amplitude
        self.randfn = randfn
        self.randargs = randargs

    def out(self, t, x):
        return self.A * self.randfn(*self.randargs)


def function(func, name="function"):
    """apply a given function to input"""
    return sl.MDLBase(None, lambda t, x, u: func(u), 1, [0], None, name=name)
