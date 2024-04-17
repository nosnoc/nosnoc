# Proposed class structure for `nosnoc`

## `nosnoc`
This is the top level namespace that all elements of nosnoc live under.
Core functionality is provided by the classes `NosnocOcp` and `NosnocIntegrator` for solving optimal control problems and initial value problems respectively.
These classes take as arguments one of the types of model and provide the functionality discussed in the class descriptions.

### `NosnocOcp`
This class provides the core optimal control functionality of nosnoc.
Constructor Arguments:
+ `model`, type: nosnoc.core.ModelBase
  + Provides model data to generate the OCP.
  + Type of model determines the kind of MPCC to generate automatically.

Properties:
+ `model`, type: nosnoc.core.ModelBase
+ `mpcc`, type: nosnoc.core.MpccBase

Methods:
+ `set(var: string, index, val)`
  + Sets the `var(index)` to the `val`.
  + more advanced setting should be available through accessing the `mpcc` property.
+ `solve()`
  + Returns: results structure for the optimal control problem.

### `NosnocIntegrator`
This class provides the core integrator functionality of nosnoc.
Constructor Arguments:
+ `model`, type: nosnoc.core.ModelBase
  + Provides model data to generate the OCP.
  + Type of model determines the kind of MPCC to generate automatically.
  
Properties:
+ `model`, type: nosnoc.core.ModelBase
+ `mpcc`, type: nosnoc.core.MpccBase

Methods
+ `step(x0, t)`
  + Returns: `x(t;x0)`

## `nosnoc.solver`

### `Options`

### `RelaxationSolver`

### `mpccsol`

## `nosnoc.solver.plugins`

## `nosnoc.core`

### `ModelBase`

### `MpccBase`

## `nosnoc.utils`

### `time_freezing`

## `nosnoc.cds`

### `Model`

### `CdsMpcc`

## `nosnoc.filippov`
Alternative name would be `nosnoc.pss`.
### `Model`

### `StewartMpcc`

### `StepMpcc`

## `nosnoc.cls`

### `Model`

### `FesdJMpcc`
