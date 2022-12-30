classdef MpccMode
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

