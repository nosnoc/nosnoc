# `nosnoc`
**nosnoc** is an open-source `MATLAB` software package for NOnSmooth Numerical Optimal Control.
The `Python` package `nosnoc_py` with similar functionality is available as well.


- Documentation is [available here](https://nosnoc.readthedocs.io/en/latest/index.html). 

- A detailed overview of the theory and methods behind nosnoc can be found in the course material of the
[Summer School on Direct Methods for Optimal Control of Nonsmooth Systems](https://www.syscop.de/teaching/ss2023/summer-school-direct-methods-optimal-control-nonsmooth-systems).


## General
**nosnoc** is a tool for numerically solving optimal control problems with nonsmooth dynamical systems that have switches and/or state jumps. It supports:

1. **Automatic discretization via the FESD method**:
   - High accuracy and correct sensitivities.
   - Classical time-stepping methods only have first-order accuracy and incorrect sensitivities, even when they appear to be differentiable.

2. **Automatic reformulations of systems with state jumps**:
   - For example, contact problems via time-freezing into Filippov systems/Piecewise Smooth Systems (PSS).
   - This enables high accuracy even for systems with state jumps.

3. **Solving nonsmooth nonlinear programs via homotopy methods**:
   - Supports numerous relaxation-based algorithms for solving mathematical programs with complementarity constraints (MPCCs). 
   - Enables the use of off-the-shelf solvers like IPOPT and SNOPT.

With **nosnoc**, one can solve optimal control problems subject to several classes of hybrid and nonsmooth systems:
- Switched systems
- Piecewise affine, piecewise smooth systems
- Filippov systems
- Systems with logical Heaviside step functions
- Relay systems
- Projected dynamical systems
- First-order Sweeping processes
- Rigid body models with impacts and friction
- Hybrid systems with hysteresis

Users may directly provide a DCS or the problem in a standard form which is automatically reformulated into a DCS.
**nosnoc** discretizes a Dynamic Complementarity System (DCS) with the FESD method and solves the resulting MPCCs. 



## Installation

**nosnoc** requires `CasADi` version >=3.5.5 and Matlab version >=R2021b, <=R2024b.

### Installation for MATLAB

1.  Install  `CasADi` and make sure it is added to your `MATLAB` path.
For `CasADi` installation instructions follow visit: https://web.casadi.org/get/

2.   Clone this repository
```
git clone --recursive https://github.com/nosnoc/nosnoc.git
```

3. Open the `nosnoc` folder in Matlab and run the `install_nosnoc` script
```
>> install_nosnoc
```

Note that `IPOPT` is shipped with `CasADi`. More information including detailed documentation can be found on its [homepage](https://coin-or.github.io/Ipopt/ )

### Installation for python

Go to the [nosnoc_py](https://github.com/nosnoc/nosnoc_py) repository for more info.

## Using nosnoc

The interface of **nosnoc** is based on the symbolic modeling framework [CasADi](https://web.casadi.org/). User inputs should be provided as `CasADi` expressions.

To get started, we recommend you look into our example libraries for:
- [MATLAB](https://github.com/nosnoc/nosnoc/tree/main/examples/matlab)
- [Python](https://github.com/nosnoc/nosnoc_py/tree/main/examples)



## Literature - theory and algorithms

### FESD
[Finite Elements with Switch Detection for Direct Optimal Control of Nonsmooth Systems](https://link.springer.com/article/10.1007/s00211-024-01412-z) \
A.Nurkanović, M. Sperl, S. Albrecht, M. Diehl \
Numerische Mathematik (2024): 1-48
```
@article{Nurkanovic2024,
  title={Finite elements with switch detection for direct optimal control of nonsmooth systems},
  author={Nurkanovi{\'c}, Armin and Sperl, Mario and Albrecht, Sebastian and Diehl, Moritz},
  journal={Numerische Mathematik},
  pages={1--48},
  year={2024},
  publisher={Springer}
}
```

[Finite Elements with Switch Detection for numerical optimal control of nonsmooth dynamical systems with set-valued heaviside step functions](https://www.sciencedirect.com/science/article/pii/S1751570X24000554) \
A. Nurkanović, A. Pozharskiy, J. Frey, M. Diehl\
Nonlinear Analysis: Hybrid Systems 54, 101518	

```
@article{Nurkanovic2024a,
  title={Finite Elements with Switch Detection for numerical optimal control of nonsmooth dynamical systems with set-valued heaviside step functions},
  author={Nurkanovi{\'c}, Armin and Pozharskiy, Anton and Frey, Jonathan and Diehl, Moritz},
  journal={Nonlinear Analysis: Hybrid Systems},
  volume={54},
  pages={101518},
  year={2024},
  publisher={Elsevier}
}
```



### Time-Freezing
[A Time-Freezing Approach for Numerical Optimal Control of Nonsmooth Differential Equations with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021.pdf) \
A. Nurkanović, T. Sartor, S. Albrecht, M. Diehl \
IEEE Control Systems Letters 2021

[The Time-Freezing Reformulation for Numerical Optimal Control of Complementarity Lagrangian Systems with State Jumps](https://www.sciencedirect.com/science/article/pii/S0005109823004594) \
A. Nurkanović, S. Albrecht, B. Brogliato, M. Diehl \
Automatica 158 (2023): 111295.

[Continuous Optimization for Control of Hybrid Systems with Hysteresis via Time-Freezing](https://cdn.syscop.de/publications/Nurkanovic2022a.pdf) \
A.Nurkanović , M. Diehl \
IEEE Control Systems Letters 2022


## Literature - Software

### nosnoc

[nosnoc: A Software Package for Numerical Optimal Control of Nonsmooth Systems](https://cdn.syscop.de/publications/Nurkanovic2022b.pdf) \
A.Nurkanović , M. Diehl \
IEEE Control Systems Letters 2022

```
@article{Nurkanovic2022,
  title={nosnoc: A software package for numerical optimal control of nonsmooth systems},
  author={Nurkanovi{\'c}, Armin and Diehl, Moritz},
  journal={IEEE Control Systems Letters},
  volume={6},
  pages={3110--3115},
  year={2022},
  publisher={IEEE}
}
```



### CasADi

[CasADi - A software framework for nonlinear optimization and optimal control](https://cdn.syscop.de/publications/Andersson2019.pdf) \
J.A.E. Andersson, J. Gillis, G. Horn, J.B. Rawlings, M. Diehl \
Mathematical Programming Computation, 2019

### IPOPT
[On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming](https://link.springer.com/article/10.1007/s10107-004-0559-y) \
A. Wächter, L. T. Biegler
Mathematical programming, 2006 

## Contact

If you have questions, remarks, or comments, you are strongly encouraged to report them by creating a new issue on this GitHub page.

Feel free to contact one of the main developers directly: 
Armin Nurkanović ([armin.nurkanovic@imtek.uni-freiburg.de](mailto:armin.nurkanovic@imtek.uni-freiburg.de)),
Anton Pozharskiy ([anton.pozharskiy@imtek.uni-freiburg.de](mailto:anton.pozharskiy@imtek.uni-freiburg.de)),
Jonathan Frey ([jonathan.frey@imtek.uni-freiburg.de](mailto:jonathan.frey@imtek.uni-freiburg.de)).

Success stories and source code contributions are very welcome.

