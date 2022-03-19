function [sol,stats,solver_initalization] = homotopy_solver(solver,model,settings,solver_initalization)
import casadi.*

%%  unfold data
unfold_struct(settings,'caller')
unfold_struct(solver_initalization,'caller')

comp_res = model.comp_res;
nabla_J_fun = model.nabla_J_fun;
s_elastic_iter = 1;

%% TODO : If complemenatiry imrpoves after an iteration this means the penalty is to low, 
% update penalty paramter and resolve the problem (this is motivated by
% Leyffer, Lopze,Nocedal Dynamic and problem ralph2) the low penalty means
% divergenc
%% Inital homotopy parameters

% if gradint_norm_penalty
%     % Leyffer paper:
%     
% else
%     % take default value, sigma_0, provided by user
% end
sigma_0_scaled = max(full(nabla_J_fun(w0,inf)));
%% data for stats

sigma_k = sigma_0;
x0 = model.x0;

complementarity_stats = [full(comp_res(w0))];
cpu_time = [];
homotopy_iterations = [];
w0_base = w0;

if store_all_homotopy_iterates
W = [w0];
end

%% homtopy loop
complementarity_iter = 1;
     ii = 0;
    while complementarity_iter > comp_tol && ii < N_homotopy && s_elastic_iter>comp_tol
        % homotopy parameter update
        if ii == 0
            sigma_k = inf;
        end
        if ii == 0
            sigma_k = sigma_0;
        end
        
        if ii>0
            sigma_k = kappa*sigma_k;
        end
        tic
        sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg,'p',sigma_k);
        cpu_time_iter = toc ;
        cpu_time = [cpu_time,cpu_time_iter];
        
        w_opt = full(sol.x);
        w0 = w_opt;
        if store_all_homotopy_iterates
            W = [W,w_opt];
        end

        % complementarity
        complementarity_iter = full(comp_res(w_opt));

        complementarity_stats = [complementarity_stats;complementarity_iter];
        % feasiblity check
        
        % primal
        
        % dual
        
        % update counter
        ii = ii+1;
        fprintf('Homotopy iteration : %d / %d, with sigma = %2.2e completed, complementarity resiudal %2.2e, CPU Iter Time:  %2.2f s .\n',ii,N_homotopy,sigma_k,complementarity_iter,cpu_time_iter);
        if mpcc_mode >= 5
            fprintf('l_inf Slack value: %2.2e  .\n',w_opt(end));
%             s_elastic_iter = w_opt(end);
        end

        if mpcc_mode >= 8       
            fprintf('rho Slack value: %2.2e  .\n',w_opt(end-1));
%             s_elastic_iter = w_opt(end);
        end

        if complementarity_iter > 1e1 && ii >= 0.75*N_homotopy
            error('The homotopy loop is diverging. Try chaning parameters of the MPCC homotopy loop.')
            break;
        end
        
    end
   fprintf('Total homotopy CPU Time : %4.4f s  ///  %2.2f min.  \n',sum(cpu_time),sum(cpu_time)/60);
%% output
stats.complementarity_stats = complementarity_stats;
stats.cpu_time = cpu_time;
stats.cpu_time_total = sum(cpu_time);
stats.w_opt = w_opt;
stats.sigma_k = sigma_k;
stats.homotopy_iterations = ii;

if store_all_homotopy_iterates
            sol.W = W;
end
end

