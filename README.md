# NOS-NOC
**NOS-NOC** is open source software package for NOnSmooth Numerical Optimal Control.


## General
It is a modular tool for numerically solving optimal control problems with piecewise smooth systems (PSS). It relies on the recently introduced Finite Elements with Switch Detection (FESD) which enables high accuracy optimal control of PSS. The time-freezing reformulation, which transforms several classes of systems with state jumps into PSS is supported as well. 
It enables the treatment of a broad class of nonsmooth systems in a unified way. The user manual can be found [here](https://github.com/nurkanovic/nosnoc/blob/main/doc/nosnoc_manual.pdf).


## Installation

**NOS-NOC** requires `CasADi` version 3.5.5.

Currently, a `MATLAB` version is avilable. Versions to come will support a `python` interface as well.
### Instalation for MATLAB


1.  Install  `CasADi` and make sure it is added to your `MATLAB` path.

     For `CasADi` installation instructions follow this [link](https://web.casadi.org/get/).
   
    
2.   Clone this repository and add it to your `MATLAB` path:

     ```
     git clone https://github.com/nurkanovic/nosnoc.git
     ```
	 

Note that `IPOPT` is shipped with `CasADi`, but more information including a detailed documentation can be found on its [homepage](https://coin-or.github.io/Ipopt/ ) 

	 
## Using NOS-NOC

The interface of **NOS-NOC** is based on the symbolic modeling framework [CasADi](https://web.casadi.org/).  
User inputs should be given as `CasADi` expressions `Function` objects.	 

Minimal code example
```matlab
import casadi.*
% Call this function to have an overview of all options.
[settings] = default_settings_fesd();  
% Highlight that the problem is treated with time-freezing reformulation and change the number of IRK stages
settings.time_freezing = 1; 
settings.n_s = 3; 
% Discretization data
model.T = 5; 
model.N_stages = 15; 
model.N_finite_elements = 3;
% Define states, controls and inital values
q = MX.sym('q',2); v = MX.sym('v',2); t = MX.sym('t');
model.x = [q;v;t];
model.x0 = [0;0.5;0;0;0];
u = MX.sym('u',2); 
model.u = u;
% Define regions of the PSS
model.S = [1; -1];
model.c = q(2); 
% Define modes of PSS 
f_1 = [v;u(1);u(2)-9.81;1];
 f_2 = [0;v(2);0;-10*(q(2))-0.211989*v(2);0];
model.F = [f_1 f_2];
% Stage and terminal costs
model.f_q = u'*u; model.f_q_T = 10*v'*v;
% Path and terminal constraints
model.g_ineq = u'*u-7^2;
model.g_terminal = q-[4;0.25];
% Solve OCP with MPCC homotopy and plot results
[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
[results,stats] = homotopy_solver(solver,model,settings,solver_initalization);
plot_result_ball(model,settings,results,stats)

````


More details can be found in the [user manual](https://github.com/nurkanovic/nosnoc/blob/main/doc/nosnoc_manual.pdf).

A `python` version is work in progress.

## Literature - theory and algortihms

### FESD
[Finite Elements with Switch Detection for Direct Optimal Control of Nonsmooth Systems](https://github.com/nurkanovic/nosnoc) \
A.Nurkanović, M. Sperl, S. Albrecht, M. Diehl \
arXiv preprint 2022

### Time - Freezing
[A Time-Freezing Approach for Numerical Optimal Control of Nonsmooth Differential Equations with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021.pdf) \
A. Nurkanović, T. Sartor, S. Albrecht, M. Diehl \
IEEE Control Systems Letters 2021

[The Time-Freezing Reformulation for Numerical Optimal Control of Complementarity Lagrangian Systems with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021a.pdf) \
A. Nurkanović, S. Albrecht, B. Brogliato, M. Diehl \
arXiv preprint 2021

[Continuous Optimization for Control of Hybrid Systems with Hysteresis via Time-Freezing](https://github.com/nurkanovic/nosnoc) \
A.Nurkanović , M. Diehl \
arXiv preprint 2022


## Literature - software

### NOS-NOC

[NOS-NOC: An Software Package for Numerical Optimal Control of Nonsmooth Systems](https://github.com/nurkanovic/nosnoc) \
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

