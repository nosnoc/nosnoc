#include<HomotopySolver.hpp>
#include<iostream>
#include<stdexcept>
#include<algorithm>

namespace nosnoc
{
HomotopySolver::HomotopySolver()
{
  casadi::Dict nlpopts;

  
    nlpopts["ipopt.print_level"] = 0;
  

  
    nlpopts["ipopt.sb"] = "yes";
  

  
    nlpopts["ipopt.max_iter"] = 500;
  

  
    nlpopts["ipopt.bound_relax_factor"] = 0;
  

  
    nlpopts["ipopt.tol"] = 1e-12;
  

  
    nlpopts["ipopt.dual_inf_tol"] = 1e-12;
  

  
    nlpopts["ipopt.compl_inf_tol"] = 1e-12;
  

  
    nlpopts["ipopt.acceptable_tol"] = 1e-06;
  

  
    nlpopts["ipopt.mu_strategy"] = "adaptive";
  

  
    nlpopts["ipopt.mu_oracle"] = "quality-function";
  

  
    nlpopts["ipopt.warm_start_init_point"] = "yes";
  

  
    nlpopts["ipopt.linear_solver"] = "mumps";
  


  
    nlpopts["print_time"] = 0;
  

  
    nlpopts["verbose"] = false;
  

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
                                               {"lam_x0", m_init_lam_w},
                                               {"lbg", m_lbg},
                                               {"ubg", m_ubg},
                                               {"lam_g0", m_init_lam_g},
                                               {"x0", w_nlp_k},
                                               {"p", m_p0}};
    
    auto res = m_nlp_solver(nlp_arg);
    w_nlp_k = res.at("x");
    
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
    ii++;
    sigma_k *= 0.1;
  } while((complementarity_iter.get_elements()[0] > 1e-09 || last_iter_failed) &&
          ii < 11 &&
          sigma_k > 1.0000000000000002e-10);

  m_w_mpcc_res = w_mpcc_res.get_elements();
  return 0;
}

std::vector<double> HomotopySolver::get(std::string var, std::vector<int> indices)
{
  std::vector<size_t> out_indices;
  
  if(var == "v_global")
  {
    std::vector<size_t> shape = {1, 1};
    std::vector<size_t> multipliers;
    for(int ii = 0-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 0-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},};
    if(indices.size() != 0)
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
  
  if(var == "T_final")
  {
    std::vector<size_t> shape = {1, 1};
    std::vector<size_t> multipliers;
    for(int ii = 0-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 0-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{1},};
    if(indices.size() != 0)
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
  
  if(var == "u")
  {
    std::vector<size_t> shape = {11, 1};
    std::vector<size_t> multipliers;
    for(int ii = 1-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 1-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{9},{55},{101},{147},{193},{239},{285},{331},{377},{423},};
    if(indices.size() != 1)
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
  
  if(var == "h")
  {
    std::vector<size_t> shape = {11, 4};
    std::vector<size_t> multipliers;
    for(int ii = 2-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 2-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{},{},{},{10},{25},{40},{},{56},{71},{86},{},{102},{117},{132},{},{148},{163},{178},{},{194},{209},{224},{},{240},{255},{270},{},{286},{301},{316},{},{332},{347},{362},{},{378},{393},{408},{},{424},{439},{454},};
    if(indices.size() != 2)
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
  
  if(var == "x")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{2, 3},{},{},{},{},{},{},{},{},{},{},{},{},{},{11, 12},{18, 19},{},{26, 27},{33, 34},{},{41, 42},{48, 49},{},{},{},{},{57, 58},{64, 65},{},{72, 73},{79, 80},{},{87, 88},{94, 95},{},{},{},{},{103, 104},{110, 111},{},{118, 119},{125, 126},{},{133, 134},{140, 141},{},{},{},{},{149, 150},{156, 157},{},{164, 165},{171, 172},{},{179, 180},{186, 187},{},{},{},{},{195, 196},{202, 203},{},{210, 211},{217, 218},{},{225, 226},{232, 233},{},{},{},{},{241, 242},{248, 249},{},{256, 257},{263, 264},{},{271, 272},{278, 279},{},{},{},{},{287, 288},{294, 295},{},{302, 303},{309, 310},{},{317, 318},{324, 325},{},{},{},{},{333, 334},{340, 341},{},{348, 349},{355, 356},{},{363, 364},{370, 371},{},{},{},{},{379, 380},{386, 387},{},{394, 395},{401, 402},{},{409, 410},{416, 417},{},{},{},{},{425, 426},{432, 433},{},{440, 441},{447, 448},{},{455, 456},{462, 463},};
    if(indices.size() != 3)
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
  
  if(var == "z")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},};
    if(indices.size() != 3)
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
  
  if(var == "lambda")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{4, 5},{},{},{},{},{},{},{},{},{},{},{},{},{},{13, 14},{20, 21},{},{28, 29},{35, 36},{},{43, 44},{50, 51},{},{},{},{},{59, 60},{66, 67},{},{74, 75},{81, 82},{},{89, 90},{96, 97},{},{},{},{},{105, 106},{112, 113},{},{120, 121},{127, 128},{},{135, 136},{142, 143},{},{},{},{},{151, 152},{158, 159},{},{166, 167},{173, 174},{},{181, 182},{188, 189},{},{},{},{},{197, 198},{204, 205},{},{212, 213},{219, 220},{},{227, 228},{234, 235},{},{},{},{},{243, 244},{250, 251},{},{258, 259},{265, 266},{},{273, 274},{280, 281},{},{},{},{},{289, 290},{296, 297},{},{304, 305},{311, 312},{},{319, 320},{326, 327},{},{},{},{},{335, 336},{342, 343},{},{350, 351},{357, 358},{},{365, 366},{372, 373},{},{},{},{},{381, 382},{388, 389},{},{396, 397},{403, 404},{},{411, 412},{418, 419},{},{},{},{},{427, 428},{434, 435},{},{442, 443},{449, 450},{},{457, 458},{464, 465},};
    if(indices.size() != 3)
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
  
  if(var == "theta")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{6, 7},{},{},{},{},{},{},{},{},{},{},{},{},{},{15, 16},{22, 23},{},{30, 31},{37, 38},{},{45, 46},{52, 53},{},{},{},{},{61, 62},{68, 69},{},{76, 77},{83, 84},{},{91, 92},{98, 99},{},{},{},{},{107, 108},{114, 115},{},{122, 123},{129, 130},{},{137, 138},{144, 145},{},{},{},{},{153, 154},{160, 161},{},{168, 169},{175, 176},{},{183, 184},{190, 191},{},{},{},{},{199, 200},{206, 207},{},{214, 215},{221, 222},{},{229, 230},{236, 237},{},{},{},{},{245, 246},{252, 253},{},{260, 261},{267, 268},{},{275, 276},{282, 283},{},{},{},{},{291, 292},{298, 299},{},{306, 307},{313, 314},{},{321, 322},{328, 329},{},{},{},{},{337, 338},{344, 345},{},{352, 353},{359, 360},{},{367, 368},{374, 375},{},{},{},{},{383, 384},{390, 391},{},{398, 399},{405, 406},{},{413, 414},{420, 421},{},{},{},{},{429, 430},{436, 437},{},{444, 445},{451, 452},{},{459, 460},{466, 467},};
    if(indices.size() != 3)
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
  
  if(var == "mu")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{8},{},{},{},{},{},{},{},{},{},{},{},{},{},{17},{24},{},{32},{39},{},{47},{54},{},{},{},{},{63},{70},{},{78},{85},{},{93},{100},{},{},{},{},{109},{116},{},{124},{131},{},{139},{146},{},{},{},{},{155},{162},{},{170},{177},{},{185},{192},{},{},{},{},{201},{208},{},{216},{223},{},{231},{238},{},{},{},{},{247},{254},{},{262},{269},{},{277},{284},{},{},{},{},{293},{300},{},{308},{315},{},{323},{330},{},{},{},{},{339},{346},{},{354},{361},{},{369},{376},{},{},{},{},{385},{392},{},{400},{407},{},{415},{422},{},{},{},{},{431},{438},{},{446},{453},{},{461},{468},};
    if(indices.size() != 3)
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
  
  if(var == "v_global")
  {
    std::vector<size_t> shape = {1, 1};
    std::vector<size_t> multipliers;
    for(int ii = 0-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 0-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},};
    if(indices.size() != 0)
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
  
  if(var == "T_final")
  {
    std::vector<size_t> shape = {1, 1};
    std::vector<size_t> multipliers;
    for(int ii = 0-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 0-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{1},};
    if(indices.size() != 0)
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
  
  if(var == "u")
  {
    std::vector<size_t> shape = {11, 1};
    std::vector<size_t> multipliers;
    for(int ii = 1-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 1-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{9},{55},{101},{147},{193},{239},{285},{331},{377},{423},};
    if(indices.size() != 1)
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
  
  if(var == "h")
  {
    std::vector<size_t> shape = {11, 4};
    std::vector<size_t> multipliers;
    for(int ii = 2-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 2-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{},{},{},{10},{25},{40},{},{56},{71},{86},{},{102},{117},{132},{},{148},{163},{178},{},{194},{209},{224},{},{240},{255},{270},{},{286},{301},{316},{},{332},{347},{362},{},{378},{393},{408},{},{424},{439},{454},};
    if(indices.size() != 2)
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
  
  if(var == "x")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{2, 3},{},{},{},{},{},{},{},{},{},{},{},{},{},{11, 12},{18, 19},{},{26, 27},{33, 34},{},{41, 42},{48, 49},{},{},{},{},{57, 58},{64, 65},{},{72, 73},{79, 80},{},{87, 88},{94, 95},{},{},{},{},{103, 104},{110, 111},{},{118, 119},{125, 126},{},{133, 134},{140, 141},{},{},{},{},{149, 150},{156, 157},{},{164, 165},{171, 172},{},{179, 180},{186, 187},{},{},{},{},{195, 196},{202, 203},{},{210, 211},{217, 218},{},{225, 226},{232, 233},{},{},{},{},{241, 242},{248, 249},{},{256, 257},{263, 264},{},{271, 272},{278, 279},{},{},{},{},{287, 288},{294, 295},{},{302, 303},{309, 310},{},{317, 318},{324, 325},{},{},{},{},{333, 334},{340, 341},{},{348, 349},{355, 356},{},{363, 364},{370, 371},{},{},{},{},{379, 380},{386, 387},{},{394, 395},{401, 402},{},{409, 410},{416, 417},{},{},{},{},{425, 426},{432, 433},{},{440, 441},{447, 448},{},{455, 456},{462, 463},};
    if(indices.size() != 3)
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
  
  if(var == "z")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},};
    if(indices.size() != 3)
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
  
  if(var == "lambda")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{4, 5},{},{},{},{},{},{},{},{},{},{},{},{},{},{13, 14},{20, 21},{},{28, 29},{35, 36},{},{43, 44},{50, 51},{},{},{},{},{59, 60},{66, 67},{},{74, 75},{81, 82},{},{89, 90},{96, 97},{},{},{},{},{105, 106},{112, 113},{},{120, 121},{127, 128},{},{135, 136},{142, 143},{},{},{},{},{151, 152},{158, 159},{},{166, 167},{173, 174},{},{181, 182},{188, 189},{},{},{},{},{197, 198},{204, 205},{},{212, 213},{219, 220},{},{227, 228},{234, 235},{},{},{},{},{243, 244},{250, 251},{},{258, 259},{265, 266},{},{273, 274},{280, 281},{},{},{},{},{289, 290},{296, 297},{},{304, 305},{311, 312},{},{319, 320},{326, 327},{},{},{},{},{335, 336},{342, 343},{},{350, 351},{357, 358},{},{365, 366},{372, 373},{},{},{},{},{381, 382},{388, 389},{},{396, 397},{403, 404},{},{411, 412},{418, 419},{},{},{},{},{427, 428},{434, 435},{},{442, 443},{449, 450},{},{457, 458},{464, 465},};
    if(indices.size() != 3)
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
  
  if(var == "theta")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{6, 7},{},{},{},{},{},{},{},{},{},{},{},{},{},{15, 16},{22, 23},{},{30, 31},{37, 38},{},{45, 46},{52, 53},{},{},{},{},{61, 62},{68, 69},{},{76, 77},{83, 84},{},{91, 92},{98, 99},{},{},{},{},{107, 108},{114, 115},{},{122, 123},{129, 130},{},{137, 138},{144, 145},{},{},{},{},{153, 154},{160, 161},{},{168, 169},{175, 176},{},{183, 184},{190, 191},{},{},{},{},{199, 200},{206, 207},{},{214, 215},{221, 222},{},{229, 230},{236, 237},{},{},{},{},{245, 246},{252, 253},{},{260, 261},{267, 268},{},{275, 276},{282, 283},{},{},{},{},{291, 292},{298, 299},{},{306, 307},{313, 314},{},{321, 322},{328, 329},{},{},{},{},{337, 338},{344, 345},{},{352, 353},{359, 360},{},{367, 368},{374, 375},{},{},{},{},{383, 384},{390, 391},{},{398, 399},{405, 406},{},{413, 414},{420, 421},{},{},{},{},{429, 430},{436, 437},{},{444, 445},{451, 452},{},{459, 460},{466, 467},};
    if(indices.size() != 3)
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
  
  if(var == "mu")
  {
    std::vector<size_t> shape = {11, 4, 3};
    std::vector<size_t> multipliers;
    for(int ii = 3-1; ii >= 0 ; ii--)
    {
      int mult = 1;
      for(int jj = 3-1; jj > ii; jj--)
      {
        mult *= shape[jj];
      }
      multipliers.push_back(mult);
    }
    std::reverse(multipliers.begin(), multipliers.end());
    std::vector<std::vector<size_t>> var_indices = {{},{},{8},{},{},{},{},{},{},{},{},{},{},{},{},{},{17},{24},{},{32},{39},{},{47},{54},{},{},{},{},{63},{70},{},{78},{85},{},{93},{100},{},{},{},{},{109},{116},{},{124},{131},{},{139},{146},{},{},{},{},{155},{162},{},{170},{177},{},{185},{192},{},{},{},{},{201},{208},{},{216},{223},{},{231},{238},{},{},{},{},{247},{254},{},{262},{269},{},{277},{284},{},{},{},{},{293},{300},{},{308},{315},{},{323},{330},{},{},{},{},{339},{346},{},{354},{361},{},{369},{376},{},{},{},{},{385},{392},{},{400},{407},{},{415},{422},{},{},{},{},{431},{438},{},{446},{453},{},{461},{468},};
    if(indices.size() != 3)
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
}