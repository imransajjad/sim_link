from sim_link.ode_solvers import *

import numpy as np
import matplotlib

matplotlib.use("Agg")  # comment this and uncomment plt.show to view test plots
import matplotlib.pyplot as plt


def test1():
    """
    since the solver expects a
    dx/dt = f(t,x,...) type function,
    define a bunch of the and test
    """

    def test_x(t, k, m, b, x0):
        x = np.array([0.0] * 2)
        wn = np.sqrt(k / m)
        gamma = b / 2 / np.sqrt(k * m)
        wd = wn * np.sqrt(1 - gamma**2)
        sigma = gamma * wn
        # s1 = np.array([-sigma + wd * 1.0j])
        # s2 = np.array([-sigma - wd * 1.0j])
        x00 = x0[0]
        v00 = x0[1]
        alpha = x00 / 2
        beta = (v00 + sigma * x00) / 2 / wd
        x[0] = 2 * (alpha * np.cos(wd * t) + beta * np.sin(wd * t)) * np.exp(-sigma * t)
        x[1] = 2 * (-alpha * wd*np.sin(wd * t) + beta *wd*np.cos(wd * t)) * np.exp(-sigma * t) + \
            2 * (alpha * np.cos(wd * t) + beta * np.sin(wd * t)) * np.exp(-sigma * t)*(-sigma)
        return x
        

    def test_f(t, x, k, m, b):
        dx = np.array([0.0] * 2)
        dx[0] = x[1]
        dx[1] = -k/m * x[0] - b/m * x[1]
        return dx
    
    def get_test_x(k,m,b,x0):
        return lambda t : test_x(t,k,m,b,x0)

    def get_test_f(k,m,b):
        return lambda t,x : test_f(t,x,k,m,b)

    test_funcs = [
        (get_test_x(2,1,0.5,np.array([2.0,0])), get_test_f(2,1,0.5), np.arange(0, 10, 0.05), "1/(s^2 + 0.5s + 2)"),
        (get_test_x(4,1,0.0,np.array([5.0,0])), get_test_f(4,1,0.0), np.arange(0, 10, 0.05), "1/(s^2 + 4)"),
        (get_test_x(8,1,-0.2,np.array([0.1,0])), get_test_f(8,1,-0.2), np.arange(0, 10, 0.05), "1/(s^2 -0.2s + 8)"),
    ] # list of x(t), dx/dt, trange, label tuples


    test_solvers = [
        (integrate, "integrate"),
        (trapz, "trapz"),
        (rungekutta4, "rk4"),
        (rungekutta4ad, "rk4ad"),
    ]

    fig, axes = plt.subplots(len(test_funcs),2, figsize=(10, 3*len(test_funcs)))
    for row, (x_actual, dx, trange, flabel) in enumerate(test_funcs):
        x0 = x_actual(trange[0])
        T0, X0 = trange, [x_actual(t) for t in trange]

        for solver, label in test_solvers:
            T1, X1 = solver(dx, trange, x0)
            axes[row,0].plot(T1, X1, label=label)
            err = [x_actual(t1)-x1 for t1,x1 in zip(T1,X1)]
            axes[row,1].plot(T1, err, label=label+f" {np.linalg.norm(err):.2e}")
        axes[row,0].plot(T0, X0, label="true")
        axes[row,0].legend()
        axes[row,1].legend()
        axes[row,0].set_title(flabel)
        axes[row,1].set_title(flabel + " error")
    
    fig.suptitle("ode_test_test1")
    fig.tight_layout()
    plt.savefig("ode_test_test1.png")
    # plt.show()


def _test2():
    x = np.arange(-50, 50, 0.01)

    e = lambda x, u: np.exp(-((x / u) ** 2)) / u
    Df = lambda t, x, u1, u2: np.sqrt(2 / np.pi * u1 * u2) * e(x, u1) * e(x, u2)

    u1 = 10.0
    U2 = np.arange(0.1, 2.51, 0.05)
    Y = 0 * U2
    for i, u2 in enumerate(U2):
        # print(u2)
        _, y = rungekutta4(Df, x, 0.0, u1, u2)
        Y[i] = y[-1]

    plt.plot(U2, np.sqrt(2 / (U2 / u1 + u1 / U2)), label="true")
    plt.plot(U2, Y, label="rk4")
    plt.legend()
    plt.savefig("ode_test_test2.png")
    # plt.show()


def test3():
    """
    since the solver expects a
    dx/dt = f(t,x,...) type function or equivalently a
    dy/dx = f(x,y,...) type function,
    ignore the second arg for integrating over domain only i.e.
    dy/dx = cos(x) for y = sin(x)
    """

    test_funcs = [
        (np.sin, lambda x, _ : np.cos(x), np.arange(-10,10, 0.01), "sin"),
        (np.exp, lambda x, _ : np.exp(x), np.arange(-10,10, 0.01), "exp"),
        (lambda x: 1/(1-np.pow(x,2)), lambda x, _ : 2*x/(1-np.pow(x,2))**2, np.arange(-10,-1.02, 0.02), "1/(1-x^2) region 1"),
        (lambda x: 1/(1-np.pow(x,2)), lambda x, _ : 2*x/(1-np.pow(x,2))**2, np.arange(-0.98,0.98, 0.02), "1/(1-x^2) region 2"),
        (lambda x: 1/(1-np.pow(x,2)), lambda x, _ : 2*x/(1-np.pow(x,2))**2, np.arange(1.02, 10, 0.02), "1/(1-x^2) region 3"),
    ] # list of y(x), dy/dx, xrange, label tuples


    test_solvers = [
        (integrate, "integrate"),
        (trapz, "trapz"),
        (rungekutta4, "rk4"),
        (rungekutta4ad, "rk4ad"),
    ]

    fig, axes = plt.subplots(len(test_funcs),2, figsize=(10, 3*len(test_funcs)))
    for row, (y_actual, dy, x, flabel) in enumerate(test_funcs):
        y0 = y_actual(x[0])
        X0, Y0 = x, y_actual(x)

        for solver, label in test_solvers:
            X1, Y1 = solver(dy, x, y0)
            axes[row,0].plot(X1, Y1, label=label)
            err = [y0-y1 for y0,y1 in zip(y_actual(X1),Y1)]
            axes[row,1].plot(X1, err, label=label+f"{np.linalg.norm(err):.2e}")
        axes[row,0].plot(X0, Y0, label="true")
        axes[row,0].legend()
        axes[row,1].legend()
        axes[row,0].set_title(flabel)
        axes[row,1].set_title(flabel + " error")

    
    fig.suptitle("ode_test_test3")
    fig.tight_layout()
    plt.savefig("ode_test_test3.png")
    # plt.show()

if __name__ == "__main__":
    test1()
    # _test2()
    test3()
