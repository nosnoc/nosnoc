function [lb,ub,g_comp] = generate_mpcc_relaxation_bounds(g_comp, settings)
    n_comp = size(g_comp, 1);
    switch settings.psi_fun_type
      case CFunctionType.SCHOLTES
        switch settings.relaxation_method
          case 'EQ'
            lb = zeros(n_comp, 1);
            ub = zeros(n_comp, 1);
          case 'INEQ'
            lb = -inf*ones(n_comp, 1);
            ub = zeros(n_comp, 1);
        end
      case CFunctionType.SCHOLTES_TWO_SIDED
        n_comp = n_comp/2;
        lb = [-inf*ones(n_comp, 1), zeros(n_comp, 1)]';
        lb = lb(:);
        ub = [zeros(n_comp, 1), inf*ones(n_comp, 1)]';
        ub = ub(:);
      case CFunctionType.FISCHER_BURMEISTER
        switch settings.relaxation_method
          case 'EQ'
            lb = zeros(n_comp, 1);
            ub = zeros(n_comp, 1);
          case 'INEQ'
            lb = -inf*ones(n_comp, 1);
            ub = zeros(n_comp, 1);
        end
      case CFunctionType.NATURAL_RESIDUAL
        switch settings.relaxation_method
          case 'EQ'
            lb = zeros(n_comp, 1);
            ub = zeros(n_comp, 1);
          case 'INEQ'
            lb = -inf*ones(n_comp, 1);
            ub = zeros(n_comp, 1);
        end
      case CFunctionType.CHEN_CHEN_KANZOW
        switch settings.relaxation_method
          case 'EQ'
            lb = zeros(n_comp, 1);
            ub = zeros(n_comp, 1);
          case 'INEQ'
            lb = -inf*ones(n_comp, 1);
            ub = zeros(n_comp, 1);
        end
      case CFunctionType.STEFFENSEN_ULBRICH
        switch settings.relaxation_method
          case 'EQ'
            lb = zeros(n_comp, 1);
            ub = zeros(n_comp, 1);
          case 'INEQ'
            lb = -inf*ones(n_comp, 1);
            ub = zeros(n_comp, 1);
        end
      case CFunctionType.STEFFENSEN_ULBRICH_POLY
        switch settings.relaxation_method
          case 'EQ'
            lb = zeros(n_comp, 1);
            ub = zeros(n_comp, 1);
          case 'INEQ'
            lb = -inf*ones(n_comp, 1);
            ub = zeros(n_comp, 1);
        end
      case CFunctionType.KANZOW_SCHWARTZ
        switch settings.relaxation_method
          case 'EQ'
            lb = zeros(n_comp, 1);
            ub = zeros(n_comp, 1);
          case 'INEQ'
            lb = -inf*ones(n_comp, 1);
            ub = zeros(n_comp, 1);
        end
      case CFunctionType.LIN_FUKUSHIMA
        n_comp = n_comp/2;
        lb = [-inf*ones(n_comp, 1), zeros(n_comp, 1)]';
        lb = lb(:);
        ub = [zeros(n_comp, 1), inf*ones(n_comp, 1)]';
        ub = ub(:);
      case CFunctionType.KADRANI
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
    end
end
