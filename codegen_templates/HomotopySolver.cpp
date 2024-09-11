#include "HomotopySolver.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

namespace nosnoc
{
HomotopySolver::HomotopySolver()
{
  casadi::Dict nlpopts;
{%for attribute in opts.opts_casadi_nlp.ipopt.keys()%}
  {% if opts.opts_casadi_nlp.ipopt[attribute] is string %}
    nlpopts["ipopt.{{attribute}}"] = "{{opts.opts_casadi_nlp.ipopt[attribute]}}";
  {% else %}
    nlpopts["ipopt.{{attribute}}"] = {{opts.opts_casadi_nlp.ipopt[attribute]}};
  {% endif %}
{% endfor %}
{% for attribute in opts.opts_casadi_nlp if attribute not in ["ipopt", "snopt", "uno", "worhp"]%}
  {% if opts.opts_casadi_nlp[attribute] is string %}
    nlpopts["{{attribute}}"] = "{{opts.opts_casadi_nlp[attribute]}}";
  {% elif opts.opts_casadi_nlp[attribute] is boolean%}
    nlpopts["{{attribute}}"] = {{opts.opts_casadi_nlp[attribute]|lower}};
  {% else %}
    nlpopts["{{attribute}}"] = {{opts.opts_casadi_nlp[attribute]}};
  {% endif %}
{% endfor %}
  m_nlp_solver = casadi::nlpsol("solver", "ipopt", "{{opts.solver_name}}_nlp.casadi", nlpopts);
  m_complementarity_function = casadi::external("comp_res", "{{opts.solver_name}}_comp.casadi");
}

uint32_t HomotopySolver::solve()
{
  m_p0.back() = 1;

  
  casadi::DM w_mpcc_res;
  casadi::DM(m_x0).get(w_mpcc_res, false, m_ind_mpcc);
  casadi::DM w_nlp_k = casadi::DM(m_x0);
  DM complementarity_iter;
  double sigma_k = {{opts.sigma_0}};
  bool last_iter_failed = false;
  uint32_t ii = 0;

  {% if opts.print_level >= 3 %}
  print_header();
  {% endif %}
  
  do
  {
    m_p0.back() = sigma_k;
    {% raw %}
    std::map<std::string, casadi::DM> nlp_arg = {{"lbx", m_lbw},
                                               {"ubx", m_ubw},
                                               {"lam_x0", m_init_lam_w},
                                               {"lbg", m_lbg},
                                               {"ubg", m_ubg},
                                               {"lam_g0", m_init_lam_g},
                                               {"x0", w_nlp_k},
                                               {"p", m_p0}};
    {% endraw %}
    auto res = m_nlp_solver(nlp_arg);
    w_nlp_k = res.at("x");
    auto stats = m_nlp_solver.stats();
    {% if opts.warm_start_duals %}
    m_init_lam_w = res.at("lam_x").get_elements();
    m_init_lam_g = res.at("lam_g").get_elements();
    {% endif %}
    w_nlp_k.get(w_mpcc_res, false, m_ind_mpcc);
    auto comp_args = {w_nlp_k, casadi::DM(m_p0)};
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

    {% if opts.print_level >= 3 %}
    print_nlp_iter(ii, sigma_k, complementarity_iter.get_elements().back(),
                   stats.at("iterations").as_dict().at("inf_pr").as_double_vector().back(),
                   stats.at("iterations").as_dict().at("inf_pr").as_double_vector().back(),
                   stats.at("iterations").as_dict().at("obj").as_double_vector().back(),
                   0.,stats.at("iter_count"), stats.at("return_status"));
    {% endif %}
    
    ii++;
    sigma_k *= {{opts.homotopy_update_slope}};
  } while((complementarity_iter.get_elements()[0] > {{opts.complementarity_tol}} || last_iter_failed) &&
          ii < {{opts.N_homotopy}} &&
          sigma_k > {{opts.sigma_N}});

  m_w_mpcc_res = w_mpcc_res.get_elements();
  return 0;
}

std::vector<double> HomotopySolver::get(std::string var, std::vector<int> indices)
{
  std::vector<size_t> out_indices;
  {% for varname in mpcc.w.variables.keys() %}
  if(var == "{{varname}}")
  {
    std::vector<size_t> shape = {{'{' ~ mpcc.w.variables[varname].ind_shape|join(', ') ~ '}'}};
    std::vector<size_t> multipliers;
    for(int ii = {{mpcc.w.variables[varname].depth}}-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = {{mpcc.w.variables[varname].depth}}-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{'{'}}{% for ii in range(mpcc.w.variables[varname].indices|length) %}{% if mpcc.w.variables[varname].indices[ii] is integer %}{{'{' ~ mpcc.w.variables[varname].indices[ii] ~ '}'}}{% else %}{{'{' ~ mpcc.w.variables[varname].indices[ii]|join(', ') ~ '}'}}{% endif %},{% endfor %}{{'}'}};
    if(indices.size() != {{mpcc.w.variables[varname].depth}})
    {
      throw std::invalid_argument("Wrong number of indices given to HomotopySolver::get");
    }
    int flat_index = 0;
    for(int ii = 0; ii < indices.size(); ii++)
    {
      flat_index += indices[ii]*multipliers[ii];
    }
    out_indices = var_indices[flat_index];
  }
  {% endfor %}
  std::vector<double> out;
  for(size_t idx : out_indices)
  {
    out.push_back(m_w_mpcc_res[idx-1]);
  }
  return out;
}


void HomotopySolver::set(std::string var, std::string field, std::vector<int> indices, std::vector<double> value)
{
  std::vector<size_t> out_indices;
  {% for varname in mpcc.w.variables.keys() %}
  if(var == "{{varname}}")
  {
    std::vector<size_t> shape = {{'{' ~ mpcc.w.variables[varname].ind_shape|join(', ') ~ '}'}};
    std::vector<size_t> multipliers;
    for(int ii = {{mpcc.w.variables[varname].depth}}-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = {{mpcc.w.variables[varname].depth}}-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{'{'}}{% for ii in range(mpcc.w.variables[varname].indices|length) %}{% if mpcc.w.variables[varname].indices[ii] is integer %}{{'{' ~ mpcc.w.variables[varname].indices[ii] ~ '}'}}{% else %}{{'{' ~ mpcc.w.variables[varname].indices[ii]|join(', ') ~ '}'}}{% endif %},{% endfor %}{{'}'}};
    if(indices.size() != {{mpcc.w.variables[varname].depth}})
    {
      throw std::invalid_argument("Wrong number of indices given to HomotopySolver::get");
    }
    int flat_index = 0;
    for(int ii = 0; ii < indices.size(); ii++)
    {
      flat_index += indices[ii]*multipliers[ii];
    }
    out_indices = var_indices[flat_index];
  }
  {% endfor %}

  if(field == "init")
  {
    int ii = 0;
    for(size_t idx : out_indices)
    {
      m_x0[idx-1] = value[ii];
      ii++;
    }
  }
  else if(field == "lb")
  {
    int ii = 0;
    for(size_t idx : out_indices)
    {
      m_lbw[idx-1] = value[ii];
      ii++;
    }
  }
  else if(field == "ub")
  {
    int ii = 0;
    for(size_t idx : out_indices)
    {
      m_ubw[idx-1] = value[ii];
      ii++;
    }
  }
  else
  {
    throw std::invalid_argument("Wrong number of indices given to HomotopySolver::get");
  }
}

void HomotopySolver::print_header()
{
  printf("\n|%-5s|%-10s|%-10s|%-10s|%-10s|%-10s|%-10s|%-10s|%-30s\n", "iter", "sigma", "compl_res", "inf_pr", "inf_du", "objective", "CPU time", "NLP iter", "status");
}

void HomotopySolver::print_nlp_iter(int ii, double sigma_k, double complementarity, double inf_pr, double inf_du, double objective, double cpu_time, int iter_count, std::string return_status)
{
  printf("|%-5d|%-10.2e|%-10.2e|%-10.2e|%-10.2e|%-10.2e|%-10.3f|%-10d|%-30s\n",
         ii, sigma_k, complementarity, inf_pr,inf_du,
         objective, cpu_time, iter_count, return_status.c_str());
}
}
