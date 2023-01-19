function [g, lbg, ubg, cost] = reformulate_complementarities(g, mpcc_mode, sigma_p, s_elastic)
% TODO make this work with a more structured problem.
% TODO mpcc_mode should be a structure with C-fun, relaxation type, etc.
    p = inputParser();
    p.FunctionName = 'reformulate_complementarities';
    
    n_g = size(g,1);
    if mpcc_mode == MpccMode.Scholtes_ineq
        g = g - sigma_p;
        ubg = zeros(n_g,1);
        lbg = -inf * ones(n_g,1);
    elseif mpcc_mode == MpccMode.Scholtes_eq
        g = g - sigma_p;
        ubg = zeros(n_g,1);
        lbg = zeros(n_g,1);
    elseif mpcc_mode == MpccMode.ell_1_penalty
        cost = sum(g);
    elseif mpcc_mode == MpccMode.elastic_ineq
        g = g - s_elastic*ones(n_g,1);
        ubg = zeros(n_g,1);
        lbg = -inf * ones(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_eq
        g = g - s_elastic*ones(n_g,1);
        ubg = zeros(n_g,1);
        lbg = zeros(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_two_sided
        g = [g-s_elastic*ones(n_g,1);g+s_elastic*ones(n_g,1)];
        ubg = [zeros(n_g,1); inf*ones(n_g,1)];
        lbg = [-inf*ones(n_g,1);  zeros(n_g,1)];
        % TODO these should be merged probably
    elseif mpcc_mode == MpccMode.elastic_ell_1_ineq
        g = g - s_elastic;
        ubg = zeros(n_g,1);
        lbg = -inf * ones(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_ell_1_eq
        g = g - s_elastic;
        ubg = zeros(n_g,1);
        lbg = zeros(n_g,1);
    elseif mpcc_mode == MpccMode.elastic_ell_1_two_sided
        g = [g-s_elastic;g+s_elastic];
        ubg = [zeros(n_g,1); inf*ones(n_g,1)];
        lbg = [-inf*ones(n_g,1);  zeros(n_g,1)];
    end
end
