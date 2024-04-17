# Proposed class structure for `nosnoc`

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
Base class for nosnoc models.

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

## `nosnoc.filippov`
Alternative name would be `nosnoc.pss`.
### `Model`

### `StewartMpcc`

### `StepMpcc`

## `nosnoc.cls`

### `Model`

### `FesdJMpcc`

### `TFMpcc`

## `nosnoc.cds`

### `Model`

### `CdsMpcc`
