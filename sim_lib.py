#!/usr/bin/python

"""A simple python library for some common functions/blocks for sim_link
most functions are written and tested to support numpy matrices"""

import numpy as np
import matplotlib

# matplotlib.use('TkAgg') # try if having trouble plotting on windows
import matplotlib.pyplot as plt
import time as time_sys

import sim_link as sl


def integrator(name="integrator"):
    def der(t, x, u):
        xdot = u
        return xdot

    def out(t, x, u):
        y = x[0:1]
        return y

    return sl.MDLBase(der, out, 1, [], 1, name=name)


def add(name="add"):
    """add two signals out=a+b"""

    def out(t, x, a, b):
        return a + b

    return sl.MDLBase(None, out, 2, [0, 1], 0, name=name)


def sub(name="sub"):
    """add two signals out=a+b"""

    def out(t, x, a, b):
        return a - b

    return sl.MDLBase(None, out, 2, [0, 1], 0, name=name)


def mul():
    """multiply two signals out=a*b"""

    # def der(t,x,a,b):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(t, x, a, b):
        return a * b

    # mul.der = der
    mul.out = out
    mul.namestring = "mul"

    mul.inargs = 2
    mul.cstates = 0
    mul.passargs = [0, 1]


def div():
    """divide two signals out=a/b"""

    # def der(t,x,a,b):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(t, x, a, b):
        return a / b

    # div.der = der
    div.out = out
    div.namestring = "div"

    div.inargs = 2
    div.cstates = 0
    div.passargs = [0, 1]


def inv():
    """inverse a signal out=1/u"""

    # def der(t,x,u):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(t, x, u):
        return 1.0 / u

    # inv.der = der
    inv.out = out
    inv.namestring = "inv"

    inv.inargs = 1
    inv.cstates = 0
    inv.passargs = [0]


def sat():
    """saturate by -1,1"""

    # def der(t,x,u):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(t, x, u):
        return u * (-1 < u) * (u < 1) + 1 * (u >= 1) - 1 * (u <= -1)

    # inv.der = der
    sat.out = out
    sat.namestring = "sat"

    sat.inargs = 1
    sat.cstates = 0
    sat.passargs = [0]


def sgn():
    """signum function"""

    # def der(t,x,u):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(t, x, u):
        return 1 * (u > 0) - 1 * (u < 0)

    # inv.der = der
    sgn.out = out
    sgn.namestring = "sgn"

    sgn.inargs = 1
    sgn.cstates = 0
    sgn.passargs = [0]


def time():
    """output time as a signal"""

    # def der(t,x,a,b):
    # 	return x[0:0]
    def out(t, x):
        return np.matrix([t])

    time.out = out
    time.namestring = "time"

    time.inargs = 0
    time.cstates = 0
    time.passargs = []


def sin():
    """sin(u)"""

    def out(t, x, u):
        return np.sin(u)

    sin.out = out
    sin.namestring = "sin"
    sin.inargs = 1
    sin.cstates = 0
    sin.passargs = [0]


def cos():
    """cos(u)"""

    def out(t, x, u):
        return np.cos(u)

    cos.out = out
    cos.namestring = "cos"
    cos.inargs = 1
    cos.cstates = 0
    cos.passargs = [0]


def exp():
    """exp(u)"""

    def out(t, x, u):
        return np.exp(u)

    exp.out = out
    exp.namestring = "exp"
    exp.inargs = 1
    exp.cstates = 0
    exp.passargs = [0]


def log():
    """log(u)"""

    def out(t, x, u):
        return np.log(u)

    log.out = out
    log.namestring = "log"
    log.inargs = 1
    log.cstates = 0
    log.passargs = [0]


class gain(object):
    """gain block out = k*u"""

    def __init__(self, k):
        self.k = k

        self.namestring = "gain (" + str(self.k) + ")"

        self.inargs = 1
        self.cstates = 0
        self.passargs = [0]

    # def der(self,t,x,u):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(self, t, x, u):
        y = self.k * u
        return y


class step(object):
    """step block out = 1*(t>ts)"""

    def __init__(self, ts):
        self.ts = ts

        self.namestring = "step (ts=" + str(self.ts) + "s)"

        self.inargs = 0
        self.cstates = 0
        self.passargs = []

    # def der(self,t,x,u):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(self, t, x):
        return np.matrix([1.0 * (t > self.ts)])


class const(object):
    """const block out = C"""

    def __init__(self, C):
        self.C = C

        self.namestring = "const (" + str(self.C) + ")"

        self.inargs = 0
        self.cstates = 0
        self.passargs = []

    # def der(self,t,x):
    # 	xdot = x[0:0]
    # 	return xdot
    def out(self, t, x):
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

    def out(self, t, x, u):
        return np.power(u, self.n)


class selector(object):
    """docstring for selector"""

    def __init__(self, switches):
        # super(selector, self).__init__()
        self.sw = switches

        self.namestring = "selector (" + str(self.sw) + ")"

        self.inargs = 1
        self.cstates = 0
        self.passargs = [0]

    def out(self, t, x, u):
        return u[self.sw]


class integrator_multi(object):
    """integrator, multi dimensional, nth order"""

    def __init__(self, ord, u_dim):
        # super(integrator, self).__init__()

        self.ord = ord

        self.cstates = ord * u_dim
        self.inargs = 1
        self.passargs = [0] if self.ord == 0 else []

        self.namestring = "integrator_multi"

    def der(self, t, x, u):
        # print(x)
        dx = np.roll(x, -1, axis=0)
        assert len(u) == u_dim
        dx[self.ord - 1 :: self.ord] = u
        return dx

    def out(self, t, x, u):
        # print(x)
        return x[0 :: self.ord]


class ss_LTI(object):
    """
    LTI system described by A,B,C,D matrices
    normal usage 	ss_LTI(A,B,C)
                                    ss_LTI(A,B,C,D)

    """

    def __init__(self, *args):

        self.A = args[0]
        self.B = args[1]
        self.C = args[2]
        if len(args) == 4:
            self.D = args[3]
            nz = np.nonzero(self.D)
            if np.size(nz) == 0:
                self.passargs = []
            else:
                self.passargs = [0]
        else:
            self.D = np.zeros((np.shape(self.C)[0], np.shape(self.B)[1]))
            self.passargs = []

        self.namestring = "ss_LTI"

        self.inargs = 1
        self.cstates = len(self.A)

    def der(self, t, x, u):
        xdot = self.A * x + self.B * u
        # print(x, xdot)
        return xdot

    def out(self, t, x, u):
        y = self.C * x + self.D * u
        return y


class fun_gen(object):
    """y = A*cos(wt+phi)+ofs"""

    defaultdata = {"A": 1.0, "omega": 2.0 * np.pi, "phi": 0.0, "bias": 0.0}

    def __init__(self, **kwargs):
        self.data = self.defaultdata.copy()
        self.data.update(kwargs)

        self.namestring = "fun_gen"

        self.inargs = 0
        self.cstates = 0
        self.passargs = []

    # def der(self,t,x):
    # 	return np.matrix([])
    def out(self, t, x):
        return np.matrix([self.data["A"] * np.cos(self.data["omega"] * t + self.data["phi"]) + self.data["bias"]])


class joystick_input(object):
    """joystick input for linux and windows"""

    def __init__(self, *args):
        self.namestring = "joystick_input"

        self.inargs = 0
        self.cstates = 0
        self.passargs = []

        from inputs import get_gamepad

        self.get_gamepad = get_gamepad

        self.memory = 0.0

    # def der(self,t,x):
    # 	return np.matrix([])

    def out(self, t, x):
        events = self.get_gamepad()
        for event in events:
            # print(event.ev_type, event.code, event.state)
            if event.code == "ABS_RZ":
                self.memory = -(1.0 * event.state - 128) / 128
        joy_out = [self.memory]

        joy_out = np.matrix(joy_out)
        return joy_out


class random(object):
    """docstring for rand"""

    def __init__(self, randfcn, *arg):
        # super(rand, self).__init__()
        self.arg = arg
        self.randfcn = randfcn
        self.namestring = "rand"

        self.inargs = 0
        self.cstates = 0
        self.passargs = []

    def out(self, t, x):
        return np.matrix(self.randfcn(*self.arg))


class function(object):
    """docstring for function"""

    def __init__(self, function, *arg):
        # super(function, self).__init__()
        self.arg = arg
        self.function = function
        self.namestring = "function" + str(function)

        self.inargs = 0
        self.cstates = 0
        self.passargs = []

    def out(self, t, x):
        return np.matrix(self.function(*self.arg))


class plot_window(object):
    """docstring for plot_window"""

    def __init__(self, cstates_i, *probe_list, **kwargs):
        super(plot_window, self).__init__()
        self.probe_list = probe_list

        P = {
            "xlegend": [str(i) for i in cstates_i],
            "ylegend": [str(p) for p in probe_list],
            "tlabel": "time (s)",
            "xlabel": "states",
            "ylabel": "outs",
            "plot_separate": True,
            "active_draw": False,
        }
        P.update(kwargs)

        self.active_draw = P["active_draw"]
        self.drawn = False
        self.start_time = 0.0

        if P["plot_separate"]:
            self.fig, (ax1, ax2) = plt.subplots(2, 1)
            self.axes = [ax1, ax2]
        else:
            print("here")
            self.fig, ax1 = plt.subplots(1, 1)
            self.axes = [ax1, ax1]
            P["xlegend"] = P["xlegend"] + P["ylegend"]
            P["ylegend"] = P["xlegend"]

            if not cstates_i:
                P["xlabel"] = ""
            if not probe_list:
                P["ylabel"] = ""
            sep = " + " if cstates_i and probe_list else ""

            P["xlabel"] = P["xlabel"] + sep + P["ylabel"]
            P["ylabel"] = P["xlabel"]

        self.states = self.axes[0].plot(*([[]] * (2 * len(cstates_i))))
        self.axes[0].legend(P["xlegend"])
        self.axes[0].set_ylabel(P["xlabel"])
        self.state_indices = list(cstates_i)

        if probe_list:
            self.outs = self.axes[1].plot(*([[]] * (2 * len(probe_list))))
            self.axes[1].legend(P["ylegend"])
            self.axes[1].set_ylabel(P["ylabel"])
            self.axes[1].set_xlabel(P["tlabel"])
        else:
            self.fig, (self.axes[0]) = plt.subplots(1, 1)
            self.axes[0].set_xlabel(P["tlabel"])

        self.Tlim = [0.0, 1.0]
        self.xlim = [-1.0, 1.0]
        self.axes[0].set_xlim(0.0, 1.0)
        self.axes[0].set_ylim(-10.0, 10.0)
        self.axes[0].grid()

        if probe_list:

            self.axes[1].set_xlim(0.0, 1.0)
            self.axes[1].set_ylim(-10.0, 10.0)
            if P["plot_separate"]:
                self.axes[1].grid()
                self.ylim = [-1.0, 1.0]
            else:
                self.ylim = self.xlim

    def return_axes(self):
        return self.fig, self.axes

    def animate(self, T, X, *args):

        if self.drawn:
            self.fig.canvas.restore_region(self.backgrounds[0])
            self.fig.canvas.restore_region(self.backgrounds[1])
        else:
            self.fig.canvas.draw()
            self.backgrounds = [self.fig.canvas.copy_from_bbox(ax.bbox) for ax in self.axes]

        self.fix_time(T)

        this_x = np.array(X)[:, self.state_indices]
        for i, trace in enumerate(self.states):
            trace.set_data(T, this_x[:, i])
            self.axes[0].draw_artist(trace)

            x_new = this_x[-1, i]
            if (x_new <= self.xlim[0]) or (x_new >= self.xlim[1]):
                self.fix_axes(this_x, 0)

        self.fig.canvas.blit(self.axes[0].bbox)

        if args and self.probe_list:
            this_y = np.array([np.array(y.probe_s(*self.probe_list)) for y in args[0]])
            for i, trace in enumerate(self.outs):
                trace.set_data(T, this_y[:, i])
                self.axes[1].draw_artist(trace)

                y_new = this_y[-1, i]
                if (y_new <= self.ylim[0]) or (y_new >= self.ylim[1]):
                    self.fix_axes(this_y, 1)

            self.fig.canvas.blit(self.axes[1].bbox)

        if not self.drawn:
            plt.pause(0.01)
            self.drawn = True
            self.start_time = time_sys.time()
        else:
            pass

    def fix_time(self, T):
        while T[-1] > self.Tlim[1]:
            if self.Tlim[1] <= 0:
                self.Tlim[1] += 1
            else:
                self.Tlim[1] = 1.5 * self.Tlim[1]
        if T[0] < self.Tlim[0]:
            self.Tlim[0] = T[0]
            self.Tlim[1] = T[-1] + 0.5
        for ax in self.axes:
            ax.set_xlim(self.Tlim[0], self.Tlim[1])

        self.fig.canvas.draw()
        self.backgrounds = [self.fig.canvas.copy_from_bbox(ax.bbox) for ax in self.axes]

        # Tmin is assumed to be non decreasing

    def fix_axes(self, this_x, ax_num):
        if ax_num == 0:
            lim = self.xlim
        elif ax_num == 1:
            lim = self.ylim
        x_min = min(np.amin(this_x[:, :]), lim[0])
        x_max = max(np.amax(this_x[:, :]), lim[1])

        lim[0] = 0.5 * x_max + 0.5 * x_min - 1.0 * (x_max - x_min)
        lim[1] = 0.5 * x_max + 0.5 * x_min + 1.0 * (x_max - x_min)
        self.axes[ax_num].set_ylim(lim[0], lim[1])

        self.fig.canvas.draw()
        self.backgrounds = [self.fig.canvas.copy_from_bbox(ax.bbox) for ax in self.axes]

    def show(self):
        # print "Total plot time, %2.3f " % (time_sys.time()-self.start_time)

        self.fig.canvas.flush_events()
        plt.show()
