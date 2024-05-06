 # Brief summary of nosnoc structure


##Settings, 
I vote for Option 1:
We have a one level Options class that has all discretization and all model-specific options (mixing them becomes in my opinion a bit messy, but we might discuss this again after extensive use of the new structure);


## model, dcs, mpcc, and solver 
After reading suggestions 1 and 2, my current vote is for a mix of the two:

Steps in using nosnoc:

+ define a nosnoc model,  e.g. model = nosnoc.model.pss(), define all expressions. (insetad of having nosnoc.core we may also have nosnoc.model.generic/nosnoc.model.base as a superclass of all concrete models)
+ internally we have a transformation of model to to a dcs, and have nosnoc.dcs.heaviside, nosnoc.dcs.stewart, nosnoc.dcs.pds, nosnoc.dcs.cls,
+ discretize the nosnoc.dcs.x and obtaind nosnoc.mpcc.x  (or nosnoc.discrete_time_problem.x) - here i am strongly in favor of doing all discretizations in a single method, i.e. a nested loop over N_FE and N_stages, and then we do the IRK equations + step eq + comp, all at the same time. 
Splitting them makes code more modular, but needs jumping back and forth to read equations, and every dcs has a slightly different fesd's, hence i believe inheriting this from a superclass will not improve readability.
++ Remark:In principle the user can directly specify a DCS by starting from dcs = nososnoc.dcs.x but in practive this wont be used i believe  - might be however useful for some experimental things. 
in parallel we discretize path and path complementarity constraints (which is trivial)
+ create a nosnoc.solver from the model (internally the transfromations, and discretization are carried out until the mpcc is created and mpccsol is called), syntax examples: 
  ++ ocp_solver = nosnoc.ocp.solver(model,problem_options,solver_options,);
  ++ integrator = nosnoc.integrator.solver(model,problem_options,solver_options);
+ solve problem via the nosnoc solver; and return structured results, use vdx internally to make a nice structure result from what mpccsol returns.  

## What's next
I suggest to first implement this for nosnoc.model.pss and stewart, and keep in parallel the old structure, hence i can during the review use both syntaxes and make a direct comparison. once this is set, we do the transfer for all other models and delete the old code from src.    


 


