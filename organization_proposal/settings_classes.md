# Settings Organization Options
This file contains an overview of the options structure for Nosnoc.

## Option 1: Single, Top level Options Class
This option is essentially just maintaining a single large pile of options which would combine every existing option minus those which particularly pertain to `MpccSolver`.
In this case we would simply have a top level `nosnoc.core.Options` with a combination of `NosnocProblemOptions` and the integrator related options from `nosnoc.solver.options`.
(Armin: I vote for this one)

## Option 2: Top Level Options Class tree structure
In this case we would have a single top level options class `nosnoc.core.Options` with several sub-options classes that form an options tree.
A possible structure of options based on the pipeline view of nosnoc:
+ `nosnoc.core.Options`
  + `dcs_options`: Contains options used for converting one of the models to a dcs including:
	+ lifting options for Heaviside step and FESD-J reformulations.
	+ FESD-J friction model/velocity lifting.
  + `time_freezing_options`: Contains options specific to the time-freezing reformulation
	+ Quadrature state for Lagrange Term.
	+ Polyhedral friction cone.
	+ Force lifting.
	+ etc.
  + `discretization_options`: Largest set of options. Includes all options used to turn a continuous time DCS into an MPCC
	+ IRK settings.
	+ Time optimality/speed of time variable use/control grid.
	+ Step equilibration options.
	+ Cross complementarity options.
	+ Constraint relaxation (terminal constraint and numerical/physical time).
	+ Initial impacts in FESD-J.
	+ Initial guesses for algebraic variables
	+ etc.
  + `integrator_options`: options used by the `nosnoc.Integrator`
	+ use previous solution as initialization
	+ real time plot
	+ store integrator step results.

These members should likely be basic classes themselves for ease of documentation/autocomplete.
