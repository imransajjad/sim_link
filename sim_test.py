import numpy as np
import ode_solvers as ode
import sim_link as sl
import sim_lib as sll
import time
import matplotlib.pyplot as plt
import scipy.integrate as scpi

def B():

    def der(t, x, b):
        xdot = np.matrix([])
        return xdot

    def out(t, x, b):
        y = 7 * b[0, 0]
        return y

    B.der = der
    B.out = out
    B.namestring = "function B"

    B.inargs = 1
    B.cstates = 0
    B.passargs = [0]


def A():

    def der(t, x, a):
        xdot = np.matrix([-4.0 * x[0, 0]])
        return xdot

    def out(t, x, a):
        y = np.matrix([x[0, 0] + a[0, 0]])
        return y

    A.der = der
    A.out = out
    A.namestring = "function A"

    A.inargs = 1
    A.cstates = 1
    A.passargs = [0]


def D():

    def der(t, x, d1, d2):
        xdot = -x
        return xdot

    def out(t, x, d1, d2):
        # print(t, x, d1, d2)
        y = np.matrix([x[0, 0] + d1[0, 0] + d2[0, 0]])
        return y

    D.der = der
    D.out = out
    D.namestring = "function D"

    D.cstates = 2
    D.inargs = 2
    D.passargs = [0, 1]


def G_old():

    def der(t, x, u):
        # x1dot = x2dot
        # x2dot = -1.0*x1 -2.0*x2 + u
        # xdot = np.array([ x[1], -1.0*x[0] -2.0*x[1] + u[0] ])
        A = np.matrix("0 1; -1 -2")
        B = np.matrix("0;1")
        xdot = A * x + B * u
        return xdot

    def out(t, x, u):
        # y = x1
        C = np.matrix("1 0")
        # print(C,type(C),x, type(x))
        y = C * x
        return y

    G_old.der = der
    G_old.out = out
    G_old.namestring = "Plant G"

    G_old.cstates = 2
    G_old.inargs = 2
    G_old.passargs = [1]


G_old()


def K_old():

    def der(t, x, x_ref, xhat):
        return np.array([])

    def out(t, x, x_ref, xhat):
        # print(x_ref[0,0], xhat)
        assert len(x_ref) == 1
        # e = np.concatenate( (np.array(x_ref),np.array([0])), axis=0)-xhat
        e = np.matrix("1 ; 0") * x_ref - xhat
        Kmat = np.matrix("1 2")
        # u = Kmat*np.reshape(e,(2,1))
        u = Kmat * e
        # u = u[0,0]
        return u
        # return np.array([u[0,0]])

    K_old.der = der
    K_old.out = out
    K_old.namestring = "Controller K"

    K_old.inargs = 1
    K_old.cstates = 0
    K_old.passargs = [0]


K_old()


def L_old():

    def der(t, x, u, y):
        return np.array([x[1], -1.0 * x[0] - 2.0 * x[1] + u[0]])

    def out(t, x, u, y):
        L = np.matrix("0.5; 1.5")
        return 1.0 * x + 0.0 * L * y

    L_old.der = der
    L_old.out = out
    L_old.namestring = "Observer L"

    L_old.inargs = 2
    L_old.cstates = 2
    L_old.passargs = [1]


L_old()


def ES():
    # external source

    def der(t, x):
        return np.matrix([])

    def out(t, x):
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

    A = np.matrix("0 1; -5 -0.01")
    B = np.matrix("0; 5")
    C = np.matrix("1 0")
    D = np.matrix("1")

    Trange = np.arange(0, 10.0, 0.001)
    sys = (sll.ss_LTI(A, B, C, D), [])
    x0 = ([0.0, 0.0], [])

    M, x0 = sl.init_MDL(sys, x0, "this")

    x0 = np.transpose(np.matrix(x0))
    # Tad,Xad = ode.rungekutta4ad(M.der, x_ref, Trange, x0, e_tol=1e-5, min_dt=1e-4 )
    T, X = ode.rungekutta4(M.der, x_ref, Trange, x0)
    Y = [M.out(t, x, x_ref) for t, x in zip(T, X)]

    print(T[-1], len(T))
    # print(Tad[-1], len(Tad)
    print(X[-1])
    print(Y[-1])
    plt.figure(1)
    plt.subplot(211)
    plt.plot(T, [np.array(x)[:, 0] for x in X])
    # plt.plot(Tad,[np.array(x)[:,0] for x in Xad])
    # plt.plot(T,[np.array(sl.go_deep(y[0]))[:,0] for y in Y] )

    plt.subplot(212)
    plt.plot(T)
    # plt.plot(Tad)

    plt.show()


def test2():
    import matplotlib.pyplot as plt

    T = np.arange(0, 20.0, 0.01)
    sys = (
        G,
        K,
        sll.gain(3.7),
        sll.add,
        sll.fun_gen(A=3.0, omega=2),
        sll.add,
        sll.gain(1.0),
        sll.function(np.random.randn, 1),
        sll.step(5.0),
        L,
        1,
        0,
    )
    x0 = ([1.0, 0.6], [], [], [], [], [], [], [0.0, 0.0], [], [])

    # sys = (sll.sub,sll.int1,sll.fun_gen(A=1.0,omega=2*np.pi*2.0,bias=1.0),sll.time)
    # x0 = ([],[0.0],[],[],[])

    # sys = (G,K,x_ref,L,0,1)
    M = sl.MDL(sys, x0, "this")
    M.print_table()

    PW = sll.plot_window([0, 1], [9, 0], [9, 1], plot_separate=False)
    # PW2 = sll.plot_window([0], [5], [2], [0], plot_separate=False, active_draw=True)
    fig, ax = PW.return_axes()
    ax[1].set_xlabel("no time?")

    x0 = np.transpose(np.matrix(M.x0))
    T, X, Y = ode.rungekutta4ad(
        M.der, T, x0, outcall=M.all_out, adaptive=False, min_dt=5e-3, e_tol=1e-4, realtime=True, plotcalls=[PW.animate]
    )

    PW.show()
    plt.plot(T)
    plt.show()

    print(T[-1])
    print(X[-1])
    # print(Y[-1])


def test4():
    import matplotlib.pyplot as plt

    T = np.arange(0, 5.0, 0.01)
    d_in = np.matrix([1.0])

    sys_1 = (G, K, sll.gain(3.7), [], L, 1, 0)
    x0_1 = ([1.0, 0.6], [], [], [], [0.0, 0.0])

    M1 = sl.MDL(sys_1, x0_1, "Model KGL1")
    # M1.print_table()

    sys_2 = (G, K, M1, [], L, 1, 0)
    x0_2 = ([1.5, 0.6], [], M1.x0, [], [0.0, 0.0], [], [])

    M2 = sl.MDL(sys_2, x0_2, "Model KGL21")
    # M2.print_table()

    sys_12 = (M1, M2, D, ES, [])
    x0_12 = (M1.x0, M2.x0, [0.0, 0.0], [], [], [])

    M = sl.MDL(sys_12, x0_12, "Model KGL1+KGL21-G")
    assert sl.verify(M)
    M.print_table()

    # M = sl.unpack_MDL(M)
    # assert sl.verify(M)
    # M.print_table()

    # M.out(0.0,x0,d_in)
    x0 = np.transpose(np.matrix(M.x0))
    T, X = ode.rungekutta4ad(M.der, d_in, T, x0)
    Y = [M.out(t, x, d_in) for t, x in zip(T, X)]

    print(T[-1])
    print(X[-1])
    print(Y[-1])
    print(Y[-1].probe_s([0], [2], [3]))
    plt.plot(T, [np.array(x)[:, 0] for x in X])
    plt.show()


def invalid_function():
    def der(t, x):
        return 0

    invalid_function.der = der

    invalid_function.cstates = 0


def test0():
    sys = (sll.gain(3.4), [])
    x0 = ([],)

    G1 = sl.MDL(sys, x0, "sys_gain1")

    M = sl.MDL((G1, sll.const(1.5)), ([], []), "sys1")
    M = sl.unpack_MDL(M)
    sys[0].k = 1.0

    # x0 = np.transpose(np.matrix(M.x0))
    t = np.arange(0, 1, 0.01)

    T, X = ode.rungekutta4ad(M.der, t, x0)
    Y = [M.out(t, x) for t, x in zip(T, X)]
    print(T[-1])
    print(X[-1])
    print(Y[-1])


def test8():
    import matplotlib.pyplot as plt

    x_ref = [1.0]

    T = np.arange(0, 10.0, 0.01)

    G = sl.MDLBase(G_old.der, G_old.out, 2, [1], [1.0, 0.6], name="G_sys")
    diff = sl.MDLBase(None, lambda a, b: a - b, 2, [0, 1], [], name="diff_sys")
    K = sl.MDLBase(K_old.der, K_old.out, 1, [0], [], name="K_sys")
    L = sl.MDLBase(L_old.der, L_old.out, 2, [1], [0.0, 0.0], name="L_sys")
    sys_cfg = (G, K, diff, [], L, 1, 0, [])
    # sys_cfg = (G,K,[],L,0,1) # try this, it'll cause an error

    M = sl.MDL(sys_cfg, "sys_model_1")
    print(M.table())
    x0 = M.get_x0()
    print(x0)

    T, X = ode.rungekutta4ad(M.der, x_ref, T, x0)
    Y = [M.out(t, x, x_ref).probe_s([0], [1]) for t, x in zip(T, X)]

    print(T[-1])
    print(X[-1])
    print(Y[-1])
    plt.plot(T, [np.array(x)[:, 0] for x in X])
    plt.show()


def test7():

    T = np.arange(-0.0, 10.0, 0.01)

    sys1 = (G, K, sll.sub(), [], L, 1, 0, sll.integrator(), 2)
    x01 = ([1.0, 0.6], [], [], [0.0, 0.0], [], [], [0.0])
    M1 = sl.MDL(sys1, "sys1")
    M1.print_table()

    sys2 = (G, K, sll.sub(), [], L, 0, 1, [])
    x02 = ([0.0, 0.0], [], [], [0.0, 0.0], [], [], [])
    M2 = sl.MDL(sys2, "sys2")
    M2.print_table()

    x_in = np.matrix("8.0")

    # sys = (sll.gain(5.2),M1,sll.sub,[],M2,0)
    # x0 = ([],[],[],[],[],[])
    # M = sl.MDL(sys, 'false algebraic loop')
    # M.print_table()
    # # M = sl.unpack_MDL(M)
    # # M.print_table()

    PW = sll.plot_window([0, 2, 4], [0], [1, 0], plot_separate=True)
    fig, axes = PW.return_axes()
    axes[1].legend(["out", "sys1-G"])

    x0 = np.transpose(np.matrix(M.x0))
    # y0 = M.all_out(T[0], x0, x_in)
    # print(y0)

    # print([ y0[i] for i in range(0,len(y0))])
    # dx0 = M.der(T[0], x0, x_in)
    # print(dx0)

    T, X, Y = ode.rungekutta4ad(
        M.der, x_in, T, x0, outcall=M.all_out, adaptive=False, realtime=False, plotcalls=[PW.animate], plottime=0.01
    )

    PW.show()

    print(T[-1])
    print(X[-1])
    print(Y[-1].probe_s([0], [1, 0], [1, 1]))


def test_GLK():
    x_ref = np.array([-1.045,0.0], ndmin=2).T
    T = np.arange(0, 10.0, 0.01)

    def G_der(t, x, u):
        A = np.matrix("0 1; -1 -2")
        B = np.matrix("0;1")
        xdot = A * x + B * u
        return xdot

    def G_out(t, x, u):
        C = np.matrix("1 0")
        y = C * x
        return y

    def K_der(t, x, e):
        return np.array([])

    def K_out(t, x, e):
        Kmat = np.matrix("1 2")
        u = Kmat * e
        return u

    def L_der(t, x, u, y):
        return np.array([x[1], -1.0 * x[0] - 2.0 * x[1] + u[0, 0]])

    def L_out(t, x, u, y):
        L = np.matrix("0.5; 1.5")
        return 1.0 * x + L * y

    G = sl.MDLBase(G_der, G_out, 1, [], np.array([1.0, 0.6], ndmin=2).T, name="G_sys")
    K = sl.MDLBase(None, K_out, 1, [0], np.array([], ndmin=2).T, name="K_sys")
    L = sl.MDLBase(L_der, L_out, 2, [1], np.array([0.0, 0.0], ndmin=2).T, name="L_sys")

    M = sl.MDL((G, K, sll.gain(10), sll.sub(), [], L, 1, 0), name="sys_model_1")
    print(M.table())
    x0 = M.get_x0()
    print("x0 =",x0)

    T, X = ode.rungekutta4ad(M.der, x_ref, T, x0)
    Y = [M.out(t, x, x_ref, probes=True) for t, x in zip(T, X)]

    print("T[-1] =",T[-1])
    print("X[-1] =",X[-1], type(X[0]))
    print("Y[-1] =",Y[-1], type(Y[0]))

    G_x = [np.array(x["G_sys"])[:,0] for x in X]
    L_x = [np.array(x["L_sys"])[:,0] for x in X]

    x_ref = [np.array(y["sys_model_1_u_0"])[:,0] for y in Y]
    G_y = [np.array(y["G_sys"])[:,0] for y in Y]
    L_y = [np.array(y["L_sys"])[:,0] for y in Y]

    fig, axes = plt.subplots(4,1)
    axes[0].plot(T, [x[0] for x in G_x], label="G_x0")
    axes[0].plot(T, [x[0] for x in L_x], label="L_x0")

    axes[1].plot(T, [x[1] for x in G_x], label="G_x1")
    axes[1].plot(T, [x[1] for x in L_x], label="L_x1")

    axes[2].plot(T, G_x, label="G_x")
    axes[2].plot(T, L_x, label="L_x")
    axes[2].plot(T, x_ref, label="x_ref")

    axes[3].plot(T, G_y, label="G_y")
    axes[3].plot(T, L_y, label="L_y")

    for ax in axes:
        ax.grid(True)
        ax.legend()

if __name__ == "__main__":
    test_GLK()
    plt.show()
