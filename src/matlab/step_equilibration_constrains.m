     
%% TODO: Make a function out of this segment
if step_equilibration && (k > 0 && (i > 0 || couple_across_stages))
            % define the switching indicator function for previous  node/finite element boundary
            if  (k > 0 && (i > 0 || couple_across_stages))
                % backward sums at stage k are equal to the forward sums at stage k-1
                sigma_lambda_B_k = sigma_lambda_F_k;
                sigma_theta_B_k = sigma_theta_F_k;
            end
            % forward sums initalized (these were computed during the defintion of the algebraic variables)
            sigma_lambda_F_k = sum_lambda_ki;
            sigma_theta_F_k = sum_theta_ki;

            if ~(i == 0 && k== 0)
                %% ???????????????? the questio is, does step eq. go across multiples shooting intervals( same is for lambda continuity)
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
                        %                     J_regularaize_delta_h  = J_regularaize_delta_h + 1*nu_k*delta_h_kq^2;
                        J_regularize_h  = J_regularize_h + tanh(nu_k/step_equilibration_sigma)*delta_h_kq^2;
                    case 2
                        J_regularize_h  = J_regularize_h + (nu_k)*delta_h_kq^2;
                    case 3
                        %                             g = {g{:}, (nu_k)*delta_h_kq};
                        %                             lbg = [lbg; 0];
                        %                             ubg = [ubg; 0];
                        M_elastic = 1e3;
                        g = {g{:}, [(nu_k)*delta_h_kq-s_elastic*M_elastic;-(nu_k)*delta_h_kq-s_elastic*M_elastic]};
                        lbg = [lbg; [-inf;-inf]];
                        ubg = [ubg; [0;0]];
                    case 4
                        %                             g = {g{:}, (nu_k)*delta_h_kq^2};
                        M_elastic = 1e4;
                        g = {g{:}, (nu_k)*delta_h_kq^2-s_elastic*M_elastic};
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
                        %                         g = {g{:}, [nu_k-nu_lift_k;nu_lift_k*delta_h_kq]};
                        %                         lbg = [lbg; 0;0];
                        %                         ubg = [ubg; 0;0];
                        M_elastic = 1e2;
                        g = {g{:}, [nu_k-nu_lift_k;nu_lift_k*delta_h_kq-s_elastic*M_elastic;-nu_lift_k*delta_h_kq-s_elastic*M_elastic]};
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

                        g = {g{:}, [tanh(nu_k/step_equilibration_sigma)-nu_lift_k;nu_lift_k*delta_h_kq^2-s_elastic]};
                        lbg = [lbg; 0;-inf];
                        ubg = [ubg; 0;0];
                    otherwise
                        error(',,,,')
                end
            end
        end