# Proposed class structure for `nosnoc`

## `nosnoc`
This is the top level namespace that all elements of nosnoc lives under.

## `nosnoc.core`

### `ModelBase`
Abstract base class for nosnoc models.
Provides shared functionality for OCPs and Simulations of the availiable models in nosnoc.

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

### `DCSBase`
Abstract Base class for DCS. Unsure if necessary as it is possible to go straight from model to MPCC.
In theory this is also impetus to move forward towards defining FESD as a discretization method for a more generic class of DCS. 
For this we essentially need the ability to have continuous and discontinuous complementarity functions.

Properties:
+ `f_x_fun`: Forward simulation function.
+ `g_path_fun`: Path constraints function
+ `g_terminal_fun`: Terminal constraints function
+ `g_path_comp_fun`: Path complementarity constraint function

Methods:
+ `generate_variables()`: Generate missing variables.
+ `generate_functions()`: Generate functions required for problem generation.

### `MpccBase`
Abstract base class for MPCCs. Subclass of `vdx.problems.Mpcc`.

Methods:
+ `create_variables()`: Creates all the variables of the Mpcc.
+ `forward_sim_constraints()`: Creates forward simulation and stagewise constraints.
+ `create_complementarities()`: Creates complementarity constraints.
+ `step_equilibration()`: Creates step equilibration constriants.

### `DiscretizationOptions`
Options class used to discretize a DCS and create an MPCC.

Properties:
+ `N_stages`:
+ `N_finite_elements`:
+ `n_s`:
+ `irk_scheme`:
+ `irk_representation`:
+ `T`:
+ `cross_comp_mode`:
+ `g_path_at_fe`:
+ `g_path_at_stage`:
+ `x_box_at_fe`:
+ `x_box_at_stg`:
+ `time_optimal_problem`:
+ `step_equiliration`, `rho_h`, `step_equilibration_sigma`:
+ `equidistant_control_grid`:
+ `use_speed_of_time_variables`, `local_speed_of_time_variable`, `stagewise_clock_constraint`, `impose_terminal_phyisical_time`:
+ `s_sot0`, `s_sot_max`, `s_sot_min`, `s_sot_nominal`, `rho_sot`:
+ `T_final_max`, `T_final_min`:
+ `lift_complementarities`:

## `nosnoc.model`

### `PSS`
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

### `AizermanPyatnitskii`
Subclass of `nosnoc.core.ModelBase`.

Properties:
+ `dims`: Structure of dimensions, contains additional fields:
  + `n_sys`: Number of independent systems.
  + `n_c_sys`: Number of switching functions.
  + `n_f_sys`: Number of smooth functions.
+ `f_x`: Smooth dynamic functions.
+ `c`: Switching functions.
+ `alpha`: Step functions.

### `CLS`
Subclass of `nosnoc.core.ModelBase`.

Properties:
+ `dims`: Structure of dimensions, contains additional fields.
  + `n_q`: dimension of generalized positions.
  + `n_contacts`: number of contacts 
+ `q`: .
+ `v`: .
+ `f_v`: .
+ `f_c`: .
+ `mu_f`: .
+ `e`:
+ `M`,
+ `J_normal`:
+ `J_tangent`:

### `CDS`
Subclass of `nosnoc.core.ModelBase`.

Properties:
+ `dims`: Structure of dimensions, contains additional fields:
  + `n_c`: Number of boundary functions.
+ `f_x`: Smooth dynamic functions.
+ `c`: Switching functions.

## `nosnoc.dcs`

### `Step`
Subclass of `nosnoc.core.DCSBase`.

Properties:
+ `model`: Underlying model.
+ `alpha`: Heaviside step function corresponding to switching functions `c`.
+ `lambda_p`, `lambda_n`: Lagrange multipliers from equality constraints in convex problem.
+ `beta`,`gamma`,`theta_step`: Lifting variables for convex multipliers for each region.
+ `f_x`: Dynamics, either generated from `F` or provided in the case of more general Aizerman–Pyatnitskii DIs.

### `Stewart`
Subclass of `nosnoc.core.DCSBase`.

Properties:
+ `model`: Underlying model.
+ `theta`: Stewart convex multipliers.
+ `lambda`: Lagrange multipliers from equality constraints in convex problem.
+ `mu`: Lagrange multipliers from inequality constraints in convex problem.

### `FESDJ`
Subclass of `nosnoc.core.DCSBase`.

Properties:
+ `model`: Underlying model.
+ `lambda_normal`:
+ `y_gap`:
+ `lambda_tangent`:
+ `gamma_d`:
+ `beta_d`:
+ `delta_d`:
+ `beta_conic`:
+ `gamma_conic`:
+ `p_vt`:
+ `n_vt`:
+ `alpha_vt`:
+ `Lambda_normal`:
+ `Y_gap`:
+ `L_vn`:
+ `Lambda_tangent`:
+ `Gamma_d`:
+ `Beta_d`:
+ `Delta_d`:
+ `Gamma`:
+ `Beta`:
+ `P_vt`:
+ `N_vt`:
+ `Alpha_vt`:
+ `z_v`:

### `GCS`
Subclass of `nosnoc.core.DCSBase`.

Properties:
+ `model`: Underlying model.

## `nosnoc.mpcc`
Namespace containing discretized MPCCs.
Each class in this namespace produces an MPCC which chan be given to `nosnoc.solver.mpccsol`.

### `Step`
Subclass of `nosnoc.core.MPCCBase`.

### `Stewart`
Subclass of `nosnoc.core.MPCCBase`.

### `FESDJ`
Subclass of `nosnoc.core.MPCCBase`.

### `GCS`
Subclass of `nosnoc.core.MPCCBase`.

## `nosnoc.solver`
A generic solver interface for MPCCs.

### `Options`
Options for the `nosnoc.solver.MpccSolver`.

### `MpccSolver`
Backend for `nosnoc.solver.mpccsol` which implements the relaxation/smoothing approaches to solving MPCCs.

### `mpccsol`
An `nlpsol` style interface for solving MPCCs.
