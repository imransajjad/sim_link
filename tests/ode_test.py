from sim_link.ode_solvers import *

import numpy as np
import matplotlib

matplotlib.use("Agg")  # comment this and uncomment plt.show to view test plots
import matplotlib.pyplot as plt


def normpdf(x, mu, sigma):
    xad = (x - mu) / sigma
    return np.exp(-(xad**2) / 2) / np.sqrt(2 * np.pi) / sigma


def fX(x):
    mu = 1.0
    sigma = 0.1
    return normpdf(x, mu, sigma)


def fY(x, X):
    mu = X
    sigma = 0.01
    return normpdf(x, mu, sigma)


def test1():
    """
    simulate the system given by m*x''(t) + b*x'(t) + k*x(t) = 0
    with initial x(0) = 2.0, x'(0) = 0.0, compare analytical solution with
    solvers output
    """
    T = np.arange(0, 10, 0.05)
    x0 = np.array([2.0, 0])

    def test_f(t, x):
        dx = np.array([0.0] * 2)

        dx[0] = x[1]
        dx[1] = -2 * x[0] - 0.5 * x[1]

        return dx

    k = 2
    m = 1
    b = 0.5
    wn = np.sqrt(k / m)
    gamma = b / 2 / np.sqrt(k * m)
    wd = wn * np.sqrt(1 - gamma**2)
    sigma = gamma * wn
    s1 = np.array([-sigma + wd * 1.0j])
    s2 = np.array([-sigma - wd * 1.0j])
    # print(k, m, b, gamma, wd, sigma, wn)

    x00 = x0[0]
    v00 = x0[1]
    alpha = x00 / 2
    beta = (v00 + sigma * x00) / 2 / wd

    X0 = [1 * np.exp(s1 * t) + 1 * np.exp(s2 * t) for t in T]
    X0 = [2 * (alpha * np.cos(wd * t) + beta * np.sin(wd * t)) * np.exp(-sigma * t) for t in T]
    T1, X1 = integrate(test_f, T, x0)
    T2, X2 = trapz(test_f, T, x0)
    T3, X3 = rungekutta4(test_f, T, x0)

    plt.plot(T, X0, label="true")
    plt.plot(T1, [x[0] for x in X1], label="int")
    plt.plot(T2, [x[0] for x in X2], label="trapz")
    plt.plot(T3, [x[0] for x in X3], label="rk4")
    plt.legend()
    plt.savefig("ode_test_test1.png")
    # plt.show()


def test2():
    x = np.arange(-50, 50, 0.01)

    e = lambda x, u: np.exp(-((x / u) ** 2)) / u
    Df = lambda t, x, u1, u2: np.sqrt(2 / np.pi * u1 * u2) * e(x, u1) * e(x, u2)

    u1 = 10.0
    U2 = np.arange(0.1, 2.51, 0.05)
    Y = 0 * U2
    for i, u2 in enumerate(U2):
        # print(u2)
        _, y = rungekutta4(Df, u1, u2, x, 0.0)
        Y[i] = y[-1]

    plt.plot(U2, np.sqrt(2 / (U2 / u1 + u1 / U2)), label="true")
    plt.plot(U2, Y, label="rk4")
    plt.legend()
    plt.savefig("ode_test_test2.png")
    # plt.show()


if __name__ == "__main__":
    test1()
    test2()
