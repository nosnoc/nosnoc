# `nosnoc`

[![Tests](https://img.shields.io/github/actions/workflow/status/nosnoc/nosnoc/ci.yml?label=tests)](https://github.com/nosnoc/nosnoc/actions)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://nosnoc.readthedocs.io/en/latest/)
[![License](https://img.shields.io/github/license/nosnoc/nosnoc)](https://github.com/nosnoc/nosnoc/blob/main/LICENSE)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b--R2024b-orange)](https://www.mathworks.com/products/matlab.html)
[![Python package](https://img.shields.io/badge/python-nosnoc__py-green)](https://github.com/nosnoc/nosnoc_py)

**nosnoc** is an open-source `MATLAB` software package for **NOnSmooth Numerical Optimal Control and Model Predictive Control** of hybrid and nonsmooth dynamical systems. 

- Documentation: [nosnoc.readthedocs.io](https://nosnoc.readthedocs.io/en/latest/index.html)
- Theory and background: [Winter School on Numerical Methods for Optimal Control of Nonsmooth Systems (with video lectures!)](https://www.syscop.de/event/winter-school-numerical-methods-optimal-control-nonsmooth-systems) and [Summer School on Direct Methods for Optimal Control of Nonsmooth Systems](https://www.syscop.de/teaching/ss2023/summer-school-direct-methods-optimal-control-nonsmooth-systems)

For a **quick start**, we recommend browsing the [`examples`](https://github.com/nosnoc/nosnoc/tree/main/examples).

A related `Python` package, [`nosnoc_py`](https://github.com/nosnoc/nosnoc_py), is available as well.

## TL;DR nosnoc is about
- optimal control and MPC for hybrid and nonsmooth systems
- FESD discretization with accurate event handling and sensitivities
- time-freezing reformulations for systems with state jumps
- real-time MPC algorithms for hybrid systems, and fast MPC via CCOpt
- extensive MATLAB example library

## General

**nosnoc** is a software framework for numerical **optimal control** and **model predictive control (MPC)** of **hybrid and nonsmooth dynamical systems**, including systems with switching, state jumps, impacts, hysteresis, and complementarity structure.

Its main capabilities include:

1. **Real-time MPC algorithms for hybrid systems**
   - hybrid real-time iterations
   - hybrid advanced-step controller
   - hybrid advanced-step real-time iteration
   - For getting started, see: [`MPC examples`](https://github.com/nosnoc/nosnoc/tree/main/examples/mpc).

   Fast MPC in **nosnoc** uses **CCOpt**, which is currently our **fastest and most robust solver**, especially suited for MPC applications.

2. **Automatic discretization via FESD (Finite Elements with Switch Detection)**
   - high accuracy and correct sensitivities
   - superior treatment of switching events compared to classical time-stepping methods

3. **Automatic reformulations of systems with state jumps**
   - for example, contact problems via **time-freezing**
   - reformulation into Filippov / piecewise smooth / complementarity-based models
   - accurate treatment of jumps and mode transitions

4. **Homotopy-based and active-set solutions of mathematical programs with complementarity constraints**
   - multiple relaxation-based algorithms for MPCCs (`mpccsol`, `CCOpt`)
   - active-set-based methods (`mpecopt`)
   - compatibility with off-the-shelf NLP solvers such as IPOPT, and SNOPT

With **nosnoc**, users can formulate and solve problems involving:

- switched systems
- rigid body models with impacts and friction (with and without time-freezing)
- piecewise affine and piecewise smooth systems
- Filippov systems
- systems with logical Heaviside step functions
- relay systems
- projected dynamical systems
- first-order sweeping processes
- hybrid systems with hysteresis

See our [`example library`](https://github.com/nosnoc/nosnoc/tree/main/examples) for a range of optimal control and MPC examples.

Users may either provide a dynamic complementarity system (DCS) directly, or formulate the problem in a standard form that is automatically reformulated into a DCS.
**nosnoc** then discretizes the resulting model with FESD and solves the resulting MPCC/NLP.

## Installation

**nosnoc** requires:

- `CasADi` version `>= 3.5.5`
- `MATLAB` version `>= R2021b, <= R2024b`

### Installation for MATLAB

1. Install [`CasADi`](https://web.casadi.org/get/) and make sure it is on your `MATLAB` path.

2. Clone this repository:

~~~bash
git clone --recursive https://github.com/nosnoc/nosnoc.git
~~~

3. Open the `nosnoc` folder in MATLAB and run:

~~~matlab
install_nosnoc
~~~

`IPOPT` is shipped with `CasADi`. More information is available on the [IPOPT homepage](https://coin-or.github.io/Ipopt/).

### Installation for Python

For Python support, see the [`nosnoc_py`](https://github.com/nosnoc/nosnoc_py) repository.

## Dependencies

### Core dependency

- [`CasADi`](https://web.casadi.org/) for symbolic modeling and derivative generation

### Recommended dependency for fast MPC and Simulation

For high-performance MPC and fastest simulation performance, **nosnoc** supports `CCOpt`:
- [`CCOpt.jl`](https://github.com/MadNLP/CCOpt.jl)

CCOpt is used by the fast MPC functionality in **nosnoc** and is currently our fastest and most robust option.

We are currently in the process of upstreaming `CCOpt` support into the next release of `CasADi`, however this is not yet complete.
As such please use the `ap/ccopt` branch of `CasADi` found in [this fork](https://github.com/apozharski/casadi).
This requires building `CasADi` from source with the CMake flags `-DWITH_CCOPT=ON -DWITH_BUILD_CCOPT=ON`, as well as running `MATLAB` with some additional environment variables.
For details please visit the README for [`libMad`](https://github.com/apozharski/libMad), the ahead-of-time compiled library containing both `MadNLP` and `CCOpt`.

## Using nosnoc

The interface of **nosnoc** is based on the symbolic modeling framework [`CasADi`](https://web.casadi.org/).  
User inputs should therefore be provided as `CasADi` expressions.

To get started, we recommend the example libraries:

- MATLAB examples: [`examples/matlab`](https://github.com/nosnoc/nosnoc/tree/main/examples)


## Citing nosnoc

If you use **nosnoc** in research, please cite the software paper:

~~~bibtex
@article{Nurkanovic2022,
  title={nosnoc: A software package for numerical optimal control of nonsmooth systems},
  author={Nurkanovi{\'c}, Armin and Diehl, Moritz},
  journal={IEEE Control Systems Letters},
  volume={6},
  pages={3110--3115},
  year={2022},
  publisher={IEEE}
}
~~~

### Recommended additional citations

Depending on which features of **nosnoc** you use, please also cite the corresponding methodological papers.

#### Real-time MPC algorithms

~~~bibtex
@Article{Nurkanovic2026a,
  Title                    = {Real-Time Algorithms for Model Predictive Control of Hybrid Dynamical Systems},
  Author                   = {Nurkanovi{\'c}, Armin and Pozharskiy, Anton and Diehl, Moritz},
  Journal                  = {arXiv preprint},
  Year                     = {2026},
  Url                      = {https://www.syscop.de/files/users/armin.nurkanovic/Nurkanovic2026a.pdf}
}
~~~

#### FESD

~~~bibtex
@article{Nurkanovic2024,
  title={Finite elements with switch detection for direct optimal control of nonsmooth systems},
  author={Nurkanovi{\'c}, Armin and Sperl, Mario and Albrecht, Sebastian and Diehl, Moritz},
  journal={Numerische Mathematik},
  pages={1--48},
  year={2024},
  publisher={Springer}
}
~~~

~~~bibtex
@article{Nurkanovic2024a,
  title={Finite Elements with Switch Detection for numerical optimal control of nonsmooth dynamical systems with set-valued heaviside step functions},
  author={Nurkanovi{\'c}, Armin and Pozharskiy, Anton and Frey, Jonathan and Diehl, Moritz},
  journal={Nonlinear Analysis: Hybrid Systems},
  volume={54},
  pages={101518},
  year={2024},
  publisher={Elsevier}
}
~~~

#### Time-freezing

~~~bibtex
@article{Nurkanovic2021,
  title={A Time-Freezing Approach for Numerical Optimal Control of Nonsmooth Differential Equations with State Jumps},
  author={Nurkanovi{\'c}, Armin and Sartor, Thomas and Albrecht, Sebastian and Diehl, Moritz},
  journal={IEEE Control Systems Letters},
  year={2021}
}
~~~

~~~bibtex
@article{Nurkanovic2023,
  title={The Time-Freezing Reformulation for Numerical Optimal Control of Complementarity Lagrangian Systems with State Jumps},
  author={Nurkanovi{\'c}, Armin and Albrecht, Sebastian and Brogliato, Bernard and Diehl, Moritz},
  journal={Automatica},
  volume={158},
  pages={111295},
  year={2023}
}
~~~

~~~bibtex
@article{Nurkanovic2022a,
  title={Continuous Optimization for Control of Hybrid Systems with Hysteresis via Time-Freezing},
  author={Nurkanovi{\'c}, Armin and Diehl, Moritz},
  journal={IEEE Control Systems Letters},
  year={2022}
}
~~~

## Related software

### CCOpt

- [CCOpt.jl](https://github.com/MadNLP/CCOpt.jl)

~~~bibtex
@Article{Pozharskiy2026,
  Title                    = {{CCO}pt: an Open-Source Solver for Large-Scale Mathematical Programs with Complementarity Constraints},
  Author                   = {Pozharskiy, Anton and Pacaud, Fran{\c{c}}ois and Diehl, Moritz and Nurkanovi{\'c}, Armin},
  Journal                  = {arXiv preprint},
  Year                     = {2026},
  Url                      = {https://www.syscop.de/files/users/armin.nurkanovic/Pozharskiy2026.pdf}
}
~~~

### CasADi

- [CasADi -- A software framework for nonlinear optimization and optimal control](https://cdn.syscop.de/publications/Andersson2019.pdf)

### IPOPT

- [On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming](https://link.springer.com/article/10.1007/s10107-004-0559-y)

## Literature

### Real-time MPC algorithms

- [Real-Time Algorithms for Model Predictive Control of Hybrid Dynamical Systems](https://www.syscop.de/files/users/armin.nurkanovic/Nurkanovic2026a.pdf)

### FESD

- [Finite Elements with Switch Detection for Direct Optimal Control of Nonsmooth Systems](https://link.springer.com/article/10.1007/s00211-024-01412-z)
- [Finite Elements with Switch Detection for numerical optimal control of nonsmooth dynamical systems with set-valued heaviside step functions](https://www.sciencedirect.com/science/article/pii/S1751570X24000554)

### Projected dynamical systems 

- [First-Order Sweeping Processes and Extended Projected Dynamical Systems: Equivalence, Time-Discretization and Numerical Optimal Control](https://publications.syscop.de/Pozharskiy2025.pdf)
- [Finite Elements with Switch Detection for Numerical Optimal Control of Projected Dynamical Systems](https://publications.syscop.de/Pozharskiy2024c.pdf)


### Time-freezing

- [A Time-Freezing Approach for Numerical Optimal Control of Nonsmooth Differential Equations with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021.pdf)
- [The Time-Freezing Reformulation for Numerical Optimal Control of Complementarity Lagrangian Systems with State Jumps](https://www.sciencedirect.com/science/article/pii/S0005109823004594)
- [Continuous Optimization for Control of Hybrid Systems with Hysteresis via Time-Freezing](https://cdn.syscop.de/publications/Nurkanovic2022a.pdf)

## Contact

Questions, remarks, bug reports, and feature requests are best submitted via a new issue in this repository.

Main developers:
- Anton Pozharskiy — [anton.pozharskiy@imtek.uni-freiburg.de](mailto:anton.pozharskiy@imtek.uni-freiburg.de)
- Armin Nurkanović — [armin.nurkanovic@imtek.uni-freiburg.de](mailto:armin.nurkanovic@imtek.uni-freiburg.de)
- Jonathan Frey — [jonathan.frey@imtek.uni-freiburg.de](mailto:jonathan.frey@imtek.uni-freiburg.de)

Success stories and source code contributions are very welcome.