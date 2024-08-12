#include<HomotopySolver.hpp>
#include<iostream>

HomotopySolver::HomotopySolver()
{
  casadi::Dict nlpopts;
  // TODO: populate all of these 
  nlpopts["ipopt.print_level"] = {{opts.opts_casadi_nlp.ipopt.print_level}};
  nlpopts["print_time"] = {{opts.opts_casadi_nlp.print_time}};
  nlpopts["ipopt.sb"] = "{{opts.opts_casadi_nlp.ipopt.sb}}";
  nlpopts["ipopt.bound_relax_factor"] = {{opts.opts_casadi_nlp.ipopt.bound_relax_factor}};
  nlpopts["ipopt.tol"] = {{opts.opts_casadi_nlp.ipopt.tol}};
  m_nlp_solver = casadi::nlpsol("solver", "ipopt", "{{opts.solver_name}}_nlp.casadi", nlpopts);
  m_complementarity_function = casadi::external("comp_res", "{{opts.solver_name}}_comp.casadi");
}

uint32_t HomotopySolver::solve(std::map<std::string, casadi::DM> arg)
{
  m_p0.back() = 1;

  
  casadi::DM w_mpcc_res;
  casadi::DM(m_x0).get(w_mpcc_res, false, m_ind_mpcc);
  casadi::DM w_nlp_k = casadi::DM(m_x0);
  DM complementarity_iter;
  double sigma_k = {{opts.sigma_0}};
  bool last_iter_failed = false;
  uint32_t ii = 0;
  do
  {
    m_p0.back() = sigma_k;
    {% raw %}
    std::map<std::string, casadi::DM> nlp_arg = {{"lbx", m_lbw},
                                               {"ubx", m_ubw},
                                               {"lbg", m_lbg},
                                               {"ubg", m_ubg},
                                               {"x0", w_nlp_k},
                                               {"p", m_p0}};
    {% endraw %}
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
    sigma_k *= {{opts.homotopy_update_slope}};
  } while((complementarity_iter.get_elements()[0] > {{opts.complementarity_tol}} || last_iter_failed) &&
          ii < {{opts.N_homotopy}} &&
          sigma_k > {{opts.sigma_N}});

  std::cout << "mpcc result w: " << w_mpcc_res.get_elements() << std::endl;
  return 0;
}
