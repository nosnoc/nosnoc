# NOSNOC
**NOSNOC** is open source software package for NOnSmooth Numerical Optimal Control.


## General
It is a modular tool for numerically solving optimal control problems with piecewise smooth systems (PSS). It relies on the recently introduced Finite Elements with Switch Detection (FESD) which enables high accuracy optimal control of PSS. The time-freezing reformulation, which transforms several classes of systems with state jumps into PSS is supported as well. 
It enables the treatment of a broad class of nonsmooth systems in a unified way. The user manual can be found [here](https://github.com/nurkanovic/nosnoc/blob/main/doc/nosnoc_manual.pdf).

## Installation

**NOSNOC** requires `CasADi` version 3.5.5.

Currently, a `MATLAB` version is avilable. Versions to come will support a `python` interface as well.
### Installation for MATLAB


1.  Install  `CasADi` and make sure it is added to your `MATLAB` path.

     For `CasADi` installation instructions follow this [link](https://web.casadi.org/get/).
   
    
2.   Clone this repository and add it to your `MATLAB` path:

     ```
     git clone https://github.com/nurkanovic/nosnoc.git
     ```
	 

Note that `IPOPT` is shipped with `CasADi`, but more information including a detailed documentation can be found on its [homepage](https://coin-or.github.io/Ipopt/ ) 

### Installation for python

A `python` version is currently under development.
	 
## Using NOSNOC

The interface of **NOSNOC** is based on the symbolic modeling framework [CasADi](https://web.casadi.org/).  
User inputs should be given as `CasADi` expressions.

Minimal code example for time-optimal problem for a car with two modes of opration.
```matlab
import casadi.*
% Call this function to have an overview of all options.
[settings] = default_settings_nosnoc();  
% Choosing the Runge - Kutta Method and number of stages
settings.irk_scheme = 'Lobatto-IIIA';
settings.n_s = 2;
% Time-settings  - Solve an time optimal control problem
settings.time_optimal_problem = 1;
% Model - define all problem functions and
% Discretization parameters
model.N_stages = 10; % number of control intervals
model.N_finite_elements = 6; % number of finite element on every control intevral (optionally a vector might be passed)
model.T = 15;    % Time horizon
% Symbolic variables and bounds
q = SX.sym('q'); v = SX.sym('v'); 
model.x = [q;v]; % add all important data to the struct model,
model.x0 = [0;0]; % inital value
% bounds on states
model.lbx = [-inf;-20];
model.ubx = [inf;20];
% control
u = SX.sym('u'); model.u = u;
model.lbu = -5; model.ubu = 5;
% Dyanmics and the regions
f_1 = [v;u]; % mode 1 - nominal
f_2 = [v;3*u]; % mode 2 - turbo
model.S = [-1;1];
model.F = [f_1 f_2];
model.c = v-10;
% Add terminal constraint
model.g_terminal = [q-200;v-0];
% Solve OCP
[results,stats,model,settings] = nosnoc_solver(model,settings);

````


More details can be found in the [user manual](https://github.com/nurkanovic/nosnoc/blob/main/doc/nosnoc_manual.pdf).



## Literature - theory and algortihms

### FESD
[Finite Elements with Switch Detection for Direct Optimal Control of Nonsmooth Systems](https://arxiv.org/abs/2205.05337) \
A.Nurkanović, M. Sperl, S. Albrecht, M. Diehl \
arXiv preprint 2022

### Time - Freezing
[A Time-Freezing Approach for Numerical Optimal Control of Nonsmooth Differential Equations with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021.pdf) \
A. Nurkanović, T. Sartor, S. Albrecht, M. Diehl \
IEEE Control Systems Letters 2021

[The Time-Freezing Reformulation for Numerical Optimal Control of Complementarity Lagrangian Systems with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021a.pdf) \
A. Nurkanović, S. Albrecht, B. Brogliato, M. Diehl \
arXiv preprint 2021

[Continuous Optimization for Control of Hybrid Systems with Hysteresis via Time-Freezing](https://arxiv.org/abs/2203.11510) \
A.Nurkanović , M. Diehl \
arXiv preprint 2022


## Literature - software

### NOSNOC

[NOSNOC: An Software Package for Numerical Optimal Control of Nonsmooth Systems](https://arxiv.org/abs/2203.11516) \
A.Nurkanović , M. Diehl \
arXiv preprint 2022



### CasADi

[CasADi - A software framework for nonlinear optimization and optimal control](https://cdn.syscop.de/publications/Andersson2019.pdf) \
J.A.E. Andersson, J. Gillis, G. Horn, J.B. Rawlings, M. Diehl \
Mathematical Programming Computation, 2019

### IPOPT
[On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming](https://link.springer.com/article/10.1007/s10107-004-0559-y) \
A. Wächter, L. T. Biegler
Mathematical programming, 2006 

## Contact

Feel free to contact the main developer directly: Armin Nurkanović, [armin.nurkanovic@imtek.uni-freiburg.de](mailto:armin.nurkanovic@imtek.uni-freiburg.de)
If you have got questions, remarks or comments, you are strongly encouraged to report them by creating a new issue on this github page. Success stories and source code contributions are very welcome.

