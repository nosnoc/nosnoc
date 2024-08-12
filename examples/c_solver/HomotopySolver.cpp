#include<HomotopySolver.hpp>
#include<iostream>

HomotopySolver::HomotopySolver()
{
  casadi::Dict nlpopts;
  // TODO: populate all of these 
  nlpopts["ipopt.print_level"] = 0;
  nlpopts["print_time"] = 0;
  nlpopts["ipopt.sb"] = "yes";
  nlpopts["ipopt.bound_relax_factor"] = 0;
  nlpopts["ipopt.tol"] = 1e-12;
  m_nlp_solver = casadi::nlpsol("solver", "ipopt", "nosnoc_solver_nlp.casadi", nlpopts);
  m_complementarity_function = casadi::external("comp_res", "nosnoc_solver_comp.casadi");
}

uint32_t HomotopySolver::solve(std::map<std::string, casadi::DM> arg)
{
  m_p0.back() = 1;

  
  casadi::DM w_mpcc_res;
  casadi::DM(m_x0).get(w_mpcc_res, false, m_ind_mpcc);
  casadi::DM w_nlp_k = casadi::DM(m_x0);
  DM complementarity_iter;
  double sigma_k = 1;
  bool last_iter_failed = false;
  uint32_t ii = 0;
  do
  {
    m_p0.back() = sigma_k;
    
    std::map<std::string, casadi::DM> nlp_arg = {{"lbx", m_lbw},
                                               {"ubx", m_ubw},
                                               {"lbg", m_lbg},
                                               {"ubg", m_ubg},
                                               {"x0", w_nlp_k},
                                               {"p", m_p0}};
    
    auto res = m_nlp_solver(nlp_arg);
    w_nlp_k = res.at("x");
    w_nlp_k.get(w_mpcc_res, false, m_ind_mpcc);
    auto p_mpcc_vec = std::vector<double>(m_p0); // TODO don't copy here for no reason
    p_mpcc_vec.pop_back();
    auto p_mpcc_dm = casadi::DM(p_mpcc_vec);
    auto comp_args = {w_mpcc_res, p_mpcc_dm};
    auto comp_fun_ret = m_complementarity_function(comp_args);
    complementarity_iter = comp_fun_ret[0];
    auto ret_status = m_nlp_solver.stats().at("return_status");
    if(std::string("Solve_Succeeded").compare(ret_status) == 0 ||
       std::string("Solved_To_Acceptable_Level").compare(ret_status) == 0||
       std::string("Search_Direction_Becomes_Too_Small").compare(ret_status) == 0)
    {
      last_iter_failed = false;
    }
    else
    {
      last_iter_failed = true;
    }
    ii++;
    sigma_k *= 0.1;
  } while((complementarity_iter.get_elements()[0] > 1e-09 || last_iter_failed) &&
          ii < 11 &&
          sigma_k > 1.0000000000000002e-10);

  std::cout << "mpcc result w: " << w_mpcc_res.get_elements() << std::endl;
  return 0;
}