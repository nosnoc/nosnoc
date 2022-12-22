# NOSNOC
**NOSNOC** is an open source software package for NOnSmooth Numerical Optimal Control.

You can use **NOSNOC**  from `MATLAB` or `python`.

## General
It is a modular tool for numerically solving nonsmooth optimal control problems with Piecewise Smooth/Filippov Systems (PSS). It supports:
1. Automatic discretization via the FESD method - high accuracy and correct sensitivities. (Note that standard time-stepping methods have only first order accuracy and wrong sensitivities even when they appear to be differentiable!)

2. Automatic reformulations of systems with state jumps (e.g. contact problems) via time-freezing into Filippov systems/PSS.
(enables high accuracy even for system with state jumps)

3. Solving the nonsmooth nonlinear programs via homotopy methods. Enables the use of off-the-shelf solvers like IPOPT.



**NOSNOC** relies on the recently introduced Finite Elements with Switch Detection (FESD) which enables high accuracy optimal control of PSS.
It enables the treatment of a broad class of nonsmooth systems in a unified way. 
The user manual can be found [here](https://github.com/nurkanovic/nosnoc/blob/main/doc/nosnoc_manual.pdf).

## Installation

**NOSNOC** requires `CasADi` version 3.5.5.

 Versions to come will support a `python` interface as well.
### Installation for MATLAB


1.  Install  `CasADi` and make sure it is added to your `MATLAB` path.

     For `CasADi` installation instructions follow this [link](https://web.casadi.org/get/).
   
    
2.   Clone this repository and add it to your `MATLAB` path:

     ```
     git clone https://github.com/nurkanovic/nosnoc.git
     ```
	 

Note that `IPOPT` is shipped with `CasADi`, but more information including a detailed documentation can be found on its [homepage](https://coin-or.github.io/Ipopt/ ) 

### Installation for python

Checkout the [submodule](https://github.com/nurkanovic/nosnoc/tree/main/external) in this repository and clone the python repository.
Afterwards,

1. Setup virtual environment:
```
virtualenv env --python=python3
```

2. Source environment:
```
source env/bin/activate
```

3. Install
```
pip install -e .
```
	 
## Using NOSNOC

The interface of **NOSNOC** is based on the symbolic modeling framework [CasADi](https://web.casadi.org/).  
User inputs should be given as `CasADi` expressions.

To get started we recommend you to check out our example library in 
[MATLAB](https://github.com/nurkanovic/nosnoc/tree/main/examples/matlab) or [python](https://github.com/FreyJo/nosnoc_py/tree/main/examples).  

In case you need help, feel free to contact us! 

More details can be found in the [user manual](https://github.com/nurkanovic/nosnoc/blob/main/doc/nosnoc_manual.pdf).



## Literature - theory and algorithms

### FESD
[Finite Elements with Switch Detection for Direct Optimal Control of Nonsmooth Systems](https://arxiv.org/abs/2205.05337) \
A.Nurkanović, M. Sperl, S. Albrecht, M. Diehl \
arXiv preprint 2022

### Time - Freezing
[A Time-Freezing Approach for Numerical Optimal Control of Nonsmooth Differential Equations with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021.pdf) \
A. Nurkanović, T. Sartor, S. Albrecht, M. Diehl \
IEEE Control Systems Letters 2021

[The Time-Freezing Reformulation for Numerical Optimal Control of Complementarity Lagrangian Systems with State Jumps](https://arxiv.org/abs/2111.06759) \
A. Nurkanović, S. Albrecht, B. Brogliato, M. Diehl \
arXiv preprint 2021

[Continuous Optimization for Control of Hybrid Systems with Hysteresis via Time-Freezing](https://cdn.syscop.de/publications/Nurkanovic2022a.pdf) \
A.Nurkanović , M. Diehl \
IEEE Control Systems Letters 2022


## Literature - software

### NOSNOC

[NOSNOC: A Software Package for Numerical Optimal Control of Nonsmooth Systems](https://cdn.syscop.de/publications/Nurkanovic2022b.pdf) \
A.Nurkanović , M. Diehl \
IEEE Control Systems Letters 2022




### CasADi

[CasADi - A software framework for nonlinear optimization and optimal control](https://cdn.syscop.de/publications/Andersson2019.pdf) \
J.A.E. Andersson, J. Gillis, G. Horn, J.B. Rawlings, M. Diehl \
Mathematical Programming Computation, 2019

### IPOPT
[On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming](https://link.springer.com/article/10.1007/s10107-004-0559-y) \
A. Wächter, L. T. Biegler
Mathematical programming, 2006 

## Contact

Feel free to contact on of the main developer directly: Armin Nurkanović, [armin.nurkanovic@imtek.uni-freiburg.de](mailto:armin.nurkanovic@imtek.uni-freiburg.de)
Jonathan Frey [jonathan.frey@imtek.uni-freiburg.de](mailto:jonathan.frey@imtek.uni-freiburg.de)
If you have got questions, remarks or comments, you are strongly encouraged to report them by creating a new issue on this github page.
Success stories and source code contributions are very welcome.

