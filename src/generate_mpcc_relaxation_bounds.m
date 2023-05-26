function [lb,ub,g_comp] = generate_mpcc_relaxation_bounds(g_comp, settings)
    n_comp = size(g_comp, 1);
    switch settings.relaxation_method
      case 'EQ'
        lb = zeros(n_comp, 1);
        ub = zeros(n_comp, 1);
      case 'INEQ'
        lb = -inf*ones(n_comp, 1);
        ub = zeros(n_comp, 1);
      case 'TWO_SIDED' % TODO this only makes sense in some cases
        n_comp = n_comp/2;
        lb = [-inf*ones(n_comp, 1), zeros(n_comp, 1)]';
        lb = lb(:);
        ub = [zeros(n_comp, 1), inf*ones(n_comp, 1)]';
        ub = ub(:);
    end
end
