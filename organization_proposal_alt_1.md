# Proposed class structure for `nosnoc`
Remaining open questions:
+ Where to put discretization options
+ Is it a good idea to implement time-freezing as a transform from one model type to another.

## `nosnoc`
This is the top level namespace that all elements of nosnoc live under.
Core functionality is provided by the classes `NosnocOcp` and `NosnocIntegrator` for solving optimal control problems and initial value problems respectively.
These classes take as arguments one of the types of model and provide the functionality discussed in the class descriptions.

### `NosnocOcp`
This class provides the core optimal control functionality of nosnoc.

Constructor Arguments:
+ `model`, type: `nosnoc.core.ModelBase`
  + Provides model data to generate the OCP.
  + Type of model determines the kind of MPCC to generate automatically.

Properties:
+ `model`, type: `nosnoc.core.ModelBase`
+ `mpcc`, type: `nosnoc.core.MpccBase`

Methods:
+ `set(var: string, index, val)`
  + Sets the `var(index)` to the `val`.
  + more advanced setting should be available through accessing the `mpcc` property.
+ `solve()`
  + Returns: results structure for the optimal control problem.

### `NosnocIntegrator`
This class provides the core integrator functionality of nosnoc.

Constructor Arguments:
+ `model`, type: `nosnoc.core.ModelBase`
  + Provides model data to generate the OCP.
  + Type of model determines the kind of MPCC to generate automatically.
  
Properties:
+ `model`, type: `nosnoc.core.ModelBase`
+ `mpcc`, type: `nosnoc.core.MpccBase`

Methods
+ `step(x0, t)`
  + Returns: `x(t;x0)`

## `nosnoc.solver`
This namespace contains a generic interface for solving MPCCs.

### `Options`
Options for the `nosnoc.solver.MpccSolver`.

### `MpccSolver`
Backend for `nosnoc.solver.mpccsol` which implements the relaxation/smoothing approaches to solving MPCCs.

### `mpccsol`
An `nlpsol` style interface for solving MPCCs.

## `nosnoc.solver.plugins`
This namespace contains plugins that handle the different options and behaviors of different NLP solvers.

## `nosnoc.core`
Contains core base classes for nosnoc models and MPCCs.

### `ModelBase`
Abstract base class for nosnoc models.

Properties:
+ `dims`: Structure of dimensions.
  + `n_x`: Dimensions of differential state.
  + `n_u`: Dimensions of controls.
  + `n_z`: Dimensions of algebraic state.
  + `n_v_global`: Dimensions of global variables.
  + `n_p_global`: Dimensions of global parameters.
  + `n_p_time_var`: Dimensions of time varying parameters.
+ `x`,`lbx`,`ubx`,`x0`: Differential states.
+ `u`,`lbu`,`ubu`,`u0`: Controls.
+ `z`,`lbz`,`ubz`,`z0`, `g_z`: Algebraic states.
+ `v_global`,`lbv_global`,`ubv_global`,`v0_global`: Global variables.
+ `p_global`, `p_global_val`: Global parameters.
+ `p_time_var`, `p_time_var_val`: Time varying parameters.
+ `g_path`, `g_path_lb`, `g_path_ub`: Path constraints.
+ `g_terminal`, `g_terminal_lb`, `g_terminal_ub`: Terminal constraints.
+ `g_comp_path`: Path complementarity constraints.
+ `f_q`: Lagrange cost term. 
+ `f_q_T`: Mayer cost term.
+ `lsq_x`,`x_ref`,`f_lsq_x`,`x_ref_val`: Running state least squares cost.
+ `lsq_u`,`u_ref`,`f_lsq_u`,`u_ref_val`: Running control least squares cost.
+ `lsq_T`,`x_ref_end`,`f_lsq_T`,`x_ref_end_val`: Terminal least squares cost.

Methods:
+ `backfill()`: Backfill missing bounds and numerical values.
+ `verify()`: Verify dimension of provided and generated variables.
+ `generate_variables()`: Generate missing variables.
+ `generate_functions()`: Generate functions required for problem generation.

### `MpccBase`
Abstract base class for MPCCs. Subclass of `vdx.problems.Mpcc`.

Methods:
+ `create_variables()`: Creates all the variables of the Mpcc.
+ `forward_sim_constraints()`: Creates forward simulation and stagewise constraints.
+ `create_complementarities()`: Creates complementarity constraints.
+ `step_equilibration()`: Creates step equilibration constriants.

## `nosnoc.filippov`
Models and problems related to piecewise-smooth dynamical systems.
An alternative name would be `nosnoc.pss`.

TODO: how should we structure the Stewart vs Heaviside Step reformulation.
It is possible that splitting into two classes may be confusing for users?
	
### `StewartModel`
Subclass of `nosnoc.core.ModelBase`.

Properties:
+ `dims`: Structure of dimensions, contains additional fields:
  + `n_sys`: Number of independent systems.
  + `n_c_sys`: Number of switching functions.
  + `n_f_sys`: Number of smooth functions.
+ `F`: Smooth dynamic functions.
+ `c`: Switching functions.
+ `S`: Sign Matrix.
+ `g_ind`: Stewart indicator functions, may be calculated from `c` and `S`.
+ `theta`: Stewart convex multipliers.
+ `lambda`: Lagrange multipliers from equality constraints in convex problem.
+ `mu`: Lagrange multipliers from inequality constraints in convex problem.

### `StewartMpcc`
Subclass of `nosnoc.core.MpccBase`.

### `StepModel`
Subclass of `nosnoc.core.ModelBase`.

Properties:
+ `dims`: Structure of dimensions, contains additional fields:
  + `n_sys`: Number of independent systems.
  + `n_c_sys`: Number of switching functions.
  + `n_f_sys`: Number of smooth functions.
+ `F`: Smooth dynamic functions.
+ `c`: Switching functions.
+ `S`: Sign Matrix.
+ `alpha`: Heaviside step function corresponding to switching functions `c`.
+ `lambda_p`, `lambda_n`: Lagrange multipliers from equality constraints in convex problem.
+ `beta`,`gamma`,`theta_step`: Lifting variables for convex multipliers for each region.
+ `f_x`: Dynamics, either generated from `F` or provided in the case of more general Aizerman–Pyatnitskii DIs.

### `StepMpcc`
Subclass of `nosnoc.core.MpccBase`.

## `nosnoc.cls`

### `Model`

### `FesdJMpcc`
Subclass of `nosnoc.core.MpccBase`.

### `time_freezing`

## `nosnoc.cds`

### `Model`
Subclass of `nosnoc.core.ModelBase`.

Properties:
+ `dims`: Structure of dimensions, contains additional fields:
  + `n_c`: Number of boundary functions.
+ `f_x`: Smooth dynamic functions.
+ `c`: Boundary functions.
+ `lambda`: Multipliers corresponding to boundary functions `c`.
+ `E`: Projection matrix for ePDS.

### `CdsMpcc`
Subclass of `nosnoc.core.MpccBase`.
