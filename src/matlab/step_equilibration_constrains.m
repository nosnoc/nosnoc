
%% TODO: Make a function out of this segment
% if step_equilibration && (k > 0 && (i > 0 || couple_across_stages))
% if step_equilibration && i < N_finite_elements(end)-1
if step_equilibration
    % define the switching indicator function for previous  node/finite element boundary
    if  (k > 0 || i > 0)
        % backward sums at current stage k are equal to the forward sums at stage previous stage (k-1)
        sigma_lambda_B_k = sigma_lambda_F_k;
        sigma_theta_B_k = sigma_theta_F_k;
    end
    % forward sums initalized (these were computed during the defintion of the algebraic variables)
    sigma_lambda_F_k = sum_Lambda_ki;
    sigma_theta_F_k = sum_Theta_ki;

    %     if ~(i == 0 && k== 0)
    %     if i>0
    if i>0 || k >0
        sigma_lambda_k = sigma_lambda_B_k+sigma_lambda_F_k;
        sigma_theta_k = sigma_theta_B_k+sigma_theta_F_k;

        pi_lambda_k = sigma_lambda_B_k.*sigma_lambda_F_k;
        pi_theta_k = sigma_theta_B_k.*sigma_theta_F_k;

        eta_k =  pi_lambda_k+ pi_theta_k;

        nu_k = 1;
        for jjj = 1:n_theta
            nu_k = nu_k*eta_k(jjj);
        end
        nu_vector = [nu_vector;nu_k];

        switch step_equilibration_mode
            case 1
                J_regularize_h  = J_regularize_h + tanh(nu_k/step_equilibration_sigma)*delta_h_ki^2;
            case 2
                J_regularize_h  = J_regularize_h + (nu_k)*delta_h_ki^2;
            case 3
                M_elastic = 1e0;
                %                 g = {g{:}, [(nu_k)*delta_h_ki-s_elastic*M_elastic;...
                %                            -(nu_k)*delta_h_ki-s_elastic*M_elastic]};
                if mpcc_mode<5
                    g = {g{:}, [tanh(nu_k/sigma_scale)*delta_h_ki-sigma*M_elastic;...
                        -tanh(nu_k/sigma_scale)*delta_h_ki-sigma*M_elastic]};
                else
                    g = {g{:}, [tanh(nu_k)*delta_h_ki-s_elastic*M_elastic;...
                        -tanh(nu_k)*delta_h_ki-s_elastic*M_elastic]};
                end
        
                lbg = [lbg; [-inf;-inf]];
                ubg = [ubg; [0;0]];
%                 g_step_eq   = [g_step_eq;[tanh(nu_k/sigma_scale)*delta_h_ki-sigma*M_elastic;...
%                         -tanh(nu_k/sigma_scale)*delta_h_ki-sigma*M_elastic]];
            case 4
                %                             g = {g{:}, (nu_k)*delta_h_ki^2};
                M_elastic = 1e0;
                g = {g{:}, tanh(nu_k/sigma_scale)*delta_h_ki^2-s_elastic*M_elastic};
                lbg = [lbg; -inf];
                ubg = [ubg; 0];
            case 5
                % same as 3 with lifted nu_k
                nu_lift_k = MX.sym('nu_lift_k', 1);
                w = {w{:}, nu_lift_k};

                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];
                w0 = [w0;1];
                lbw = [lbw; -inf];
                ubw = [ubw; inf];
                %                         g = {g{:}, [nu_k-nu_lift_k;nu_lift_k*delta_h_ki]};
                %                         lbg = [lbg; 0;0];
                %                         ubg = [ubg; 0;0];
                M_elastic = 1e2;
                g = {g{:}, [nu_k-nu_lift_k;nu_lift_k*delta_h_ki-s_elastic*M_elastic;-nu_lift_k*delta_h_ki-s_elastic*M_elastic]};
                lbg = [lbg; 0;-inf;-inf];
                ubg = [ubg; 0;0;0];
            case 6
                % same as 4 with lifted nu_k
                nu_lift_k = MX.sym('nu_lift_k', 1);
                w = {w{:}, nu_lift_k};
                w0 = [w0;1];
                ind_total  = [ind_total,ind_total(end)+1:ind_total(end)+1];

                lbw = [lbw; -inf];
                ubw = [ubw; inf];

                g = {g{:}, [tanh(nu_k/step_equilibration_sigma)-nu_lift_k;nu_lift_k*delta_h_ki^2-s_elastic]};
                lbg = [lbg; 0;-inf];
                ubg = [ubg; 0;0];
            otherwise
                error(',,,,')
        end
    end
end