#include <casadi/casadi.hpp>
#include <stdint.h>
#include <vector>
using namespace casadi;

class HomotopySolver{
 public:
  HomotopySolver();
  uint32_t solve(std::map<std::string, DM> arg);
 private:
  std::vector<double> m_lbw = {{'{' ~ nlp_lbw|join(', ') ~ '}' }};
  std::vector<double> m_ubw = {{'{' ~ nlp_ubw|join(', ') ~ '}' }};
  std::vector<double> m_lbg = {{'{' ~ nlp_lbg|join(', ') ~ '}' }};
  std::vector<double> m_ubg = {{'{' ~ nlp_ubg|join(', ') ~ '}' }};
  std::vector<double> m_p0 = {{'{' ~ nlp_p0|join(', ') ~ '}' }};
  std::vector<double> m_x0 = {{'{' ~ nlp_x0|join(', ') ~ '}' }};
  const std::vector<int> m_ind_mpcc = {{'{' ~ ind_mpcc|join(', ') ~ '}' }};
  Function m_nlp_solver;
  Function m_complementarity_function;
};
