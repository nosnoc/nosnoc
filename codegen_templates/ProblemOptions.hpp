#include<vector>
namespace nosnoc
{
namespace problem_opts
{
constexpr size_t N_stages = {{N_stages}};
const std::vector<size_t> N_finite_elements = {% if N_finite_elements is integer%}{{'{' ~ [0, N_finite_elements] ~ '}'}}{% else %}{{'{' ~ ([0] + N_finite_elements)|join(',' ) ~ '}'}}{% endif %};
constexpr size_t n_s = {{n_s}};
}
}
