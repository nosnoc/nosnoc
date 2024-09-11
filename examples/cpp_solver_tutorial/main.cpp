#include"HomotopySolver.hpp"
#include"ProblemOptions.hpp"
#include<iostream>

int main(){
  nosnoc::HomotopySolver solver;

  std::map<std::string, casadi::DM> nlp_arg = {};
  solver.set("x", "lb", {0,0,nosnoc::problem_opts::n_s}, {0, 5});
  solver.set("x", "ub", {0,0,nosnoc::problem_opts::n_s}, {0, 5});
  solver.set("x", "init", {0,0,nosnoc::problem_opts::n_s}, {0, 5});
  solver.solve(nlp_arg);

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
