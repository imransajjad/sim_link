# sim_link

## Introduction

This is my project to replicate some Simulink (R) functionality in python.

sim_link.py is a module that can be imported, it has extensive documentation in its docstring. It provides an interface to write system blocks, and connect them together to build models. These models can further be connected into higher level models and can also be unpacked into the base system blocks.

## Installation

You can use pip with the github repo

```
pip install git+https://github.com/imransajjad/sim_link.git@main
```

or clone the repo and install
```
git clone https://github.com/imransajjad/sim_link.git
pip install sim_link
```
and/or pass the ```-e``` flag to install in editable mode
```
pip install -e sim_link
```

## Basic workflow:

Step 1: Write functions representing models of dynamic systems and then create
MDL-like objects with the supplied ```MDLBase``` class, for example:

```python
import sim_link
import numpy as np

def L_der(t, x, u, y):
    return np.array([x[1], -1.0 * x[0] - 2.0 * x[1] + u[0, 0]])

def L_out(t, x, u, y):
    L = np.matrix("0.5; 1.5")
    return 1.0 * x + L * y

L = sim_link.MDLBase(L_der, L_out, 2, [1], np.array([0.0, 0.0], ndmin=2).T, name="Observer L")

# 2 is number of inputs (u and y)
# [1] are the indices of inputs that are required to compute output (y)
# np.array([0.0, 0.0]) is the initial state vector of this block
# name="Observer L" is a human readable name to be used to identify the block
```


Step 2: Connect such systems together and evaluate

When providing a tuple to the system think of it as unpacking an expression for the output as a function of the inputs:

```
G(K(diff(x_ref,L(K(...),G(...)))))  => (G,K,diff,[],L,1,0)
```

elements in this tuple can only be:
```
function objects 	<=>	system blocks
empty lists []		<=>	inputs
ints e.g 0, 1		<=> references to other signals
```
The order of elements here is important!!!

For example, here is test code for the system:
```
          __     ------         ------
x_ref -->|+-|-->|  K   |------>|  G   |-------->
          --     ------   |     ------   |
          ^               V              |
          |             ------           |
           ------------|  L   |<---------
                        ------
```
If this is taken to be the usual plant, observer controller system,
the signals on the wires should be evaluated in the order G, L, K, which is accomplished by the following code

```python
import sim_link
import sim_link.sim_lib as sim_lib
import sim_link.ode_solvers as ode
import numpy as np
import matplotlib.pyplot as plt

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

G = sim_link.MDLBase(G_der, G_out, 1, [], np.array([1.0, 0.6], ndmin=2).T, name="G_sys")
K = sim_link.MDLBase(None, K_out, 1, [0], np.array([], ndmin=2).T, name="K_sys")
L = sim_link.MDLBase(L_der, L_out, 2, [1], np.array([0.0, 0.0], ndmin=2).T, name="L_sys")

M = sim_link.MDL((G, K, sim_lib.Gain(10), sim_lib.sub(), [], L, 1, 0), name="sys_model_1")
# M = sim_link.MDL((G, K, sim_lib.gain(10), sim_lib.sub(), [], L, 0, 1), name="sys_model_1") # try this, it'll cause an algebraic loop error

print(M.table())

# M is now a model like object and M.der can be passed to an ODE solver
x_ref = np.array([-1.045, 0.0], ndmin=2).T
T = np.arange(0, 10.0, 0.01)
T, X = ode.rungekutta4ad(M.der, x_ref, T, M.get_x0())
Y = [M.out(t, x, x_ref, probes=True) for t, x in zip(T, X)]

print("T[-1] =", T[-1])
print("X[-1] =", X[-1], type(X[0]))
print("Y[-1] =", Y[-1], type(Y[0]))

G_x = [np.array(x["G_sys"])[:, 0] for x in X]
L_x = [np.array(x["L_sys"])[:, 0] for x in X]

x_ref = [np.array(y["sys_model_1_u_0"])[:, 0] for y in Y]
G_y = [np.array(y["G_sys"])[:, 0] for y in Y]
L_y = [np.array(y["L_sys"])[:, 0] for y in Y]

fig, axes = plt.subplots(4, 1)
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

plt.show()
```

There are other examples in the ```sim_link/sim_test.py``` file.


## Basic ideas

1. t,x are global (in a sense)
    * so construct a ```der(t,x,*args)``` and ```out(t,x,*args)``` by connecting individual functions, and 'extract' hidden states from subsystems into the global state vector

2. Evaluation order has to be correct for output
    * no algebraic loops are allowed,
    * MDL builds an execution table which der and out then use
    * order of der execution doesn't really matter as long as out has been evaluated
    * passargs are necessary to eliminate false positives in algebraic loops

3. This is similar to a Moore State machine
    * The results of all systems are stored in a register which is then used
    * by subsequent systems, in the proper order

4. Need so solve the differential equation formed by der and out:
    * ```xdot = der(t,x,*args)```
    * ```y = out(t,x,*args)```
    * You can use someone else's solver, write your own, or use the ode_solvers
    module

5. The MDL function is also similar to G, L etc:

    It has ```MDL(), MDL.out(t,x), MDL.der(t,x)``` methods, which are used in the same way. In addition, it has members that contain references to the original functions, execution order and valid flags, that together form the execution table.

6. Multiple instances of the same function (block) are allowed:

    Since the state is not stored inside the functions themselves, 'two or more' of the same system will act as two separate systems and their states and results will be stored separately, for example if ```A(s) = 1/s``` then
        
    ```A(A(in)) => (A,A,[in]) => 1/s^2```
    
    If you want to use the result of a function twice, use references (ints)

7. Since MDL objects also act as dynamical system functions:

    They can also be connected together
    ```python
    sys_1 = (G,K,[],L,1,0)
    M1,x0_1 = sl.MDL(sys_1, x0_1, "Model KGL1")

    sys_2 = (G,K,[],L,1,0)
    M2,x0_2 = sl.MDL(sys_2, x0_2, "Model KGL2")

    sys_12 = (M1,M2,[])
    M,x0 = sl.MDL(sys_12 ,x0_12, "Model KGL1*KGL2")
    ```
    is a valid operation

8. Minimal Dependencies:
    sim_link.py only has copy as dependencies and as such, leaves issues
    of dimension, valid numerical and mathematical operation to the user. It does
    however detect algebraic loops. It has been tested with functions that used
    numpy. Run the tests in sim_test.py to get an impression of the workflow.

## Extras

An ode_solvers.py module that provides solvers using trapezoidal, Runge-Kutta methods is also provided.
