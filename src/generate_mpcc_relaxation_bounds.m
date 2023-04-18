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
        error('Not Implemented');
    end
end
