function [lb,ub,g_comp] = generate_mpcc_relaxation_bounds(g_comp, relaxation_type)
    import nosnoc.reg_homotopy.*
    n_comp = size(g_comp, 1);
    switch MpccMethod(relaxation_type)
      case MpccMethod.SCHOLTES_EQ
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.SCHOLTES_INEQ
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.FISCHER_BURMEISTER_EQ
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.FISCHER_BURMEISTER_INEQ
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.NATURAL_RESIDUAL_EQ
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.NATURAL_RESIDUAL_INEQ
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.CHEN_CHEN_KANZOW_EQ
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.CHEN_CHEN_KANZOW_INEQ
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.STEFFENSEN_ULBRICH_EQ
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.STEFFENSEN_ULBRICH_INEQ
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.STEFFENSEN_ULBRICH_POLY_EQ
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.STEFFENSEN_ULBRICH_POLY_INEQ
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.KANZOW_SCHWARTZ_EQ
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.KANZOW_SCHWARTZ_INEQ
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case MpccMethod.LIN_FUKUSHIMA
        n_comp = n_comp/2;
        lb = [-inf*ones(n_comp, 1), zeros(n_comp, 1)]';
        lb = lb(:);
        ub = [zeros(n_comp, 1), inf*ones(n_comp, 1)]';
        ub = ub(:);
      case MpccMethod.KADRANI
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
    end
end
