#include <casadi/casadi.hpp>
#include <stdint.h>
#include <vector>
using namespace casadi;

class HomotopySolver{
 public:
  HomotopySolver();
  uint32_t solve(std::map<std::string, DM> arg);
 private:
  std::vector<double> m_lbw = {0, 0};
  std::vector<double> m_ubw = {inf, inf};
  std::vector<double> m_lbg = {-inf};
  std::vector<double> m_ubg = {0};
  std::vector<double> m_p0 = {0, 1};
  std::vector<double> m_x0 = {0, 0};
  const std::vector<int> m_ind_mpcc = {0, 1};
  Function m_nlp_solver;
  Function m_complementarity_function;
};