# NOS-NOC
**NOS-NOC** is open source software package for NOnSmooth Numerical Optimal Control.
It is a modular tool for numerically solving optimal control problems with piecewise smooth systems (PSS). It relies on the recently introduced Finite Elements with Switch Detection (FESD) which enables high accuracy optimal control of PSS. The time-freezing reformulation, which transforms several classes of systems with state jumps into PSS is supported as well. 
Hence, it enables the treatment a broad class of nonsmooth systems in a unified way. The algorithms and reformulation yield nonsmooth nonlinear programs (NLP), which can be solved with techniques of continuous optimization in a homotopy procedure, without the use of integer variables.

## Installation

**NOS-NOC** requires CasADi versions 3.5.5.

Currently, a `MATLAB` version is avilable. Versions to come will support a `python` interface as well.
### Instalation for MATLAB


x = ca.MX.sym('x',<insert state dimension>) # states

1.  Install  `CasADi` and make sure it is added to your `MATLAB` path:

     For installation instructions follow this [link](https://web.casadi.org/get/)
   
    
2.   Clone this repository and add it to your `MATLAB` path:

     ```
     git clone https://github.com/nurkanovic/nosnoc.git
     ```
	 
	 
## Using NOS-NOC

The interface of **NOS-NOC** is based on the symbolic modeling framework [CasADi](https://web.casadi.org/).  
User inputs should be given as `CasADi` expressions `Function` objects.	 

minimal code example
```matlab
import casadi.* 
x = var
	[settings] = default_settings_fesd();  	
	settings.N_homotopy = 15;
	model = michael_schumacher_problem();
		[solver,solver_initalization, model,settings] = create_nlp_fesd(model,settings);
	[results,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization);
	plot_results(results,stats)
````

## Literature - theory and algortihms

## FESD
[Finite Elements with Switch Detection for Direct Optimal Control of Nonsmooth Systems](https://github.com/nurkanovic/nosnoc) \
A.Nurkanovic, M. Sperl, S. Albrecht, M. Diehl \
arXiv preprint 2022

## Time - Freezing
[A Time-Freezing Approach for Numerical Optimal Control of Nonsmooth Differential Equations with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021.pdf) \
A. Nurkanovic, T. Sartor, S. Albrecht, M. Diehl \
IEEE Control Systems Letters 2021

[The Time-Freezing Reformulation for Numerical Optimal Control of Complementarity Lagrangian Systems with State Jumps](https://cdn.syscop.de/publications/Nurkanovic2021a.pdf) \
A. Nurkanovic, S. Albrecht, B. Brogliato, M. Diehl \
arXiv preprint 2021

[Continuous Optimization for Control of Hybrid Systems with Hysteresis via Time-Freezing](https://github.com/nurkanovic/nosnoc) \
A.Nurkanovic , M. Diehl \
arXiv preprint 2022


## Software

### NOS-NOC

[NOS-NOC: An Software Package for Numerical Optimal Control of Nonsmooth Systems](https://github.com/nurkanovic/nosnoc) \
A.Nurkanovic , M. Diehl \
arXiv preprint 2022



### CasADi

[CasADi - A software framework for nonlinear optimization and optimal control](https://cdn.syscop.de/publications/Andersson2019.pdf) \
J.A.E. Andersson, J. Gillis, G. Horn, J.B. Rawlings, M. Diehl \
Mathematical Programming Computation, 2019

### IPOPT
[On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming](https://link.springer.com/article/10.1007/s10107-004-0559-y) \
A. Wächter, L. T. Biegler
Mathematical programming, 2006
`IPOPT` is shipped with `CasADi`, but more information including a detailed documentation can be found on its [Homepage](https://coin-or.github.io/Ipopt/ ) \

## Contact

Feel free to contact the main author directly: Armin Nurkanovic, [armin.nurkanovic@imtek.uni-freiburg.de](mailto:armin.nurkanovic@imtek.uni-freiburg.de)

If you have got questions, remarks or comments, you are strongly encouraged to report them by creating a new issue on this github page.

Bug reports, success stories, source code contributions are very welcome.

## Acknowledgments

This research was supported by the DFG via Research Unit FOR 2401 and project 424107692 and by the EU via ELO-X 953348. 

