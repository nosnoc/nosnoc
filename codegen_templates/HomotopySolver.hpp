#include <casadi/casadi.hpp>
#include <stdint.h>
#include <vector>
using namespace casadi;

namespace nosnoc
{
class HomotopySolver{
 public:
  HomotopySolver();
  uint32_t solve(std::map<std::string, DM> arg);
  std::vector<double> get(std::string var, std::vector<int> indices);
  void set(std::string var, std::string field, std::vector<int> indices, std::vector<double> value);
 private:
  std::vector<double> m_lbw = {{'{' ~ nlp_lbw|join(', ') ~ '}' }};
  std::vector<double> m_ubw = {{'{' ~ nlp_ubw|join(', ') ~ '}' }};
  std::vector<double> m_init_lam_w = {{'{' ~ nlp_init_lam_w|join(', ') ~ '}' }};
  std::vector<double> m_lbg = {{'{' ~ nlp_lbg|join(', ') ~ '}' }};
  std::vector<double> m_ubg = {{'{' ~ nlp_ubg|join(', ') ~ '}' }};
  std::vector<double> m_init_lam_g = {{'{' ~ nlp_init_lam_g|join(', ') ~ '}' }};
  std::vector<double> m_p0 = {{'{' ~ nlp_p0|join(', ') ~ '}' }};
  std::vector<double> m_x0 = {{'{' ~ nlp_x0|join(', ') ~ '}' }};
  std::vector<double> m_w_mpcc_res;
  const std::vector<int> m_ind_mpcc = {{'{' ~ ind_mpcc|join(', ') ~ '}' }};
  Function m_nlp_solver;
  Function m_complementarity_function;
};
}
