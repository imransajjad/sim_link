# sim_link

This is my project to replicate some Simulink (R) functionality in python.

sim_link.py is a module that can be imported, it has extensive documentation in
its docstring. It provides an interface to write system blocks, and connect 
them together to build models. These models can further be connected into 
higher level models and can also be unpacked into the base system blocks.

sim_link.py only has copy, traceback as dependencies and as such, leaves issues
of dimension, valid numerical and mathematical operation to the user. It does
however detect algebraic loops. It has been tested with functions that used
numpy. Run the tests in sim_test.py to get an impression of the workflow.

It uses the provided ode_solvers.py module that provides solvers using
trapezoidal, Runge-Kutta methods.