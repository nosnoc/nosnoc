classdef MpccMode
% Brief MPCC Wiki
% There are several possible MPCC Solution strategies avilable, by setting mpcc_mode to :
% 'direct' - treat complementarity conditions directly in the NLP, the bilinear term is tread as an inequality constraint.
% 'Scholtes_eq' - Smooth the complementarity conditions, Scholtes' smoothing.
% 'Scholtes_ineq' - Relax the complementarity conditions, Scholtes' relaxation.
% 'ell_1_penalty' - \ell_1 penalty, penalize the sum of all bilinear terms in the objective
% 'elastic_ineq' - \ell_infty elastic mode, upper bound all bilinear term with a positive slack, and penalize the slack in the objective.
% 'elastic_eq' - \ell_infty elastic mode, equate all bilinear term to a positive slack, and penalize the slack in the objective.
% 'elastic_two_sided' - \ell_infty, same as 'elastic_ineq' but two sided.
% 'elastic_ell_1_ineq' - \ell_1, elastic mode but penalize ell_1 norm of complementarities
% 'elastic_ell_1_eq' - \ell_1, elastic mode but penalize ell_1 norm of complementarities
% 'elastic_ell_1_two_sided' - \ell_1, elastic mode but penalize ell_1 norm of complementarities

    enumeration
        Scholtes_ineq
        Scholtes_eq
        ell_1_penalty
        elastic_ineq
        elastic_eq
        elastic_two_sided
        elastic_ell_1_ineq
        elastic_ell_1_eq
        elastic_ell_1_two_sided
    end

    properties(Constant)
        elastic_ell_1 = ["elastic_ell_1_ineq"
                         "elastic_ell_1_eq"
                         "elastic_ell_1_two_sided"];
        elastic = ["elastic_ineq"
                   "elastic_eq"
                   "elastic_two_sided"];
    end
    
end

