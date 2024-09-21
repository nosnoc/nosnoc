#include"HomotopySolver.hpp" // Include the homotopy solver header for the homotopy solver class.
#include"ProblemOptions.hpp" // Optionally also include the problem options header which contains discretization info.
#include<iostream>

int main(){
  // Instantiation of the homotopy solver.
  // All data regarding discretization, indices, bounds, etc. are created at generation time (see HomotopySolver.* for how this is done)
  nosnoc::HomotopySolver solver;

  // Now we can update the bounds and initialization of the primal variables by name and index.
  // Recall in this case that for a nosnoc problem using a RADAU-IIA discretization the initial state value
  // is at the index (0,0,n_s). Here we modify the initial state of the system from [0,0] to [0,5].
  solver.set("x", "lb", {0,0,nosnoc::problem_opts::n_s}, {0, 0});
  solver.set("x", "ub", {0,0,nosnoc::problem_opts::n_s}, {0, 0});
  solver.set("x", "init", {0,0,nosnoc::problem_opts::n_s}, {0, 0});

  // Now we can call the solve() method on the solver to run the homotopy procedure.
  // This populates the internal solution vector in the HomotopySolver object.
  // It also returns zero iff the solution was successful, and non-zero otherwise.
  int ret_val = solver.solve();

  // We can now get the results using the `HomotopySolver::get()` method which takes a variable and an index vector.
  // Here we first print the control outputs and then the primal states at the end of each control stage.
  for(int ii = 1; ii <= nosnoc::problem_opts::N_stages; ii++)
  {
    std::cout << "u(" << ii << "): " << solver.get("u", {ii}) << std::endl;
  }
  
  for(int ii = 0; ii <= nosnoc::problem_opts::N_stages; ii++)
  {
    int fe_end = nosnoc::problem_opts::N_finite_elements[ii];
    std::cout << "x(" << ii << "): " << solver.get("x", {ii, fe_end, nosnoc::problem_opts::n_s}) << std::endl;
  }

  // Now we resolve with a new set of initial conditions and display the output
  std::cout << "----------------------------------------------------------------------------------------";
  int fe_end = nosnoc::problem_opts::N_finite_elements[1];
  auto init_point = solver.get("x", {1, fe_end, nosnoc::problem_opts::n_s});
  solver.set("x", "lb", {0,0,nosnoc::problem_opts::n_s}, init_point);
  solver.set("x", "ub", {0,0,nosnoc::problem_opts::n_s}, init_point);
  solver.set("x", "init", {0,0,nosnoc::problem_opts::n_s}, init_point);
  int new_ret_val = solver.solve();

  for(int ii = 1; ii <= nosnoc::problem_opts::N_stages; ii++)
  {
    std::cout << "u(" << ii << "): " << solver.get("u", {ii}) << std::endl;
  }
  
  for(int ii = 0; ii <= nosnoc::problem_opts::N_stages; ii++)
  {
    int fe_end = nosnoc::problem_opts::N_finite_elements[ii];
    std::cout << "x(" << ii << "): " << solver.get("x", {ii, fe_end, nosnoc::problem_opts::n_s}) << std::endl;
  }
}
