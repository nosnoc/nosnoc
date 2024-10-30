classdef TestClsToCls < matlab.unittest.TestCase
    methods(Test)
        function test_basic(tc)
            import casadi.*

            model = tc.simplest_model();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Cls(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.lambda_normal, [1,1]);
            tc.verifySize(dcs.y_gap, [1,1]);
            tc.verifySize(dcs.Lambda_normal, [1,1]);
            tc.verifySize(dcs.Y_gap, [1,1]);
            
            dcs.generate_equations(opts);
            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '[v, (-9.81+lambda_normal)]');
            tc.verifyEqual(str(dcs.g_alg), '(y_gap-q)');
            tc.verifyEqual(str(dcs.g_impulse), '[((v_post_impact-v_pre_impact)-Lambda_normal), (Y_gap-q), ((P_vn-N_vn)-(v_post_impact+v_pre_impact))]');
        end

        function test_carts(tc)
            import casadi.*

            model = tc.carts();
            opts = nosnoc.Options();
            opts.preprocess();
            model.verify_and_backfill(opts);

            dcs = nosnoc.dcs.Cls(model);
            dcs.generate_variables(opts);

            tc.verifySize(dcs.lambda_normal, [5,1]);
            tc.verifySize(dcs.y_gap, [5,1]);
            tc.verifySize(dcs.lambda_tangent, [5,1]);
            tc.verifySize(dcs.gamma, [5,1]);
            tc.verifySize(dcs.beta, [5,1]);
            tc.verifySize(dcs.p_vt, [5,1]);
            tc.verifySize(dcs.n_vt, [5,1]);
            tc.verifySize(dcs.Lambda_normal, [5,1]);
            tc.verifySize(dcs.Y_gap, [5,1]);
            tc.verifySize(dcs.Lambda_tangent, [5,1]);
            tc.verifySize(dcs.Gamma, [5,1]);
            tc.verifySize(dcs.Beta, [5,1]);
            tc.verifySize(dcs.P_vt, [5,1]);
            tc.verifySize(dcs.N_vt, [5,1]);
            
            dcs.generate_equations(opts);
            % Using string comparison for now.
            % In future it would be nice to be able to check equality
            % via a walk down the expression graph but that is future work.
            tc.verifyEqual(str(dcs.f_x), '@1=-0.707107, @2=-9.81, @3=0.707107, [v_0, v_1, v_2, v_3, v_4, v_5, ((@1*lambda_normal_0)+lambda_tangent_2), ((@2+lambda_normal_2)+((@1*lambda_tangent_0)+(@1*lambda_tangent_1))), ((u+((@3*lambda_normal_0)+(@1*lambda_normal_1)))+lambda_tangent_3), ((@2+lambda_normal_3)+(@3*lambda_tangent_0)), ((@3*lambda_normal_1)+lambda_tangent_4), ((@2+lambda_normal_4)+(@3*lambda_tangent_1))]');
            tc.verifyEqual(str(dcs.g_alg), '@1=1, @2=0.707107, @3=-0.707107, @4=2, @5=0.1, @6=0.2, [(y_gap_0-(((q_2-q_0)-@1)-@1)), (y_gap_1-(((q_4-q_2)-@1)-@1)), (y_gap_2-(q_1-@1)), (y_gap_3-(q_3-@1)), (y_gap_4-(q_5-@1)), (((@2*v_1)+(@3*v_3))-((@4*gamma_0)*lambda_tangent_0)), (beta_0-(sq((@5*lambda_normal_0))-sq(fabs(lambda_tangent_0)))), (((@3*v_1)+(@2*v_3))-(p_vt_0-n_vt_0)), (((@2*v_1)+(@3*v_5))-((@4*gamma_1)*lambda_tangent_1)), (beta_1-(sq((@5*lambda_normal_1))-sq(fabs(lambda_tangent_1)))), (((@3*v_1)+(@2*v_5))-(p_vt_1-n_vt_1)), (-(v_0+((@4*gamma_2)*lambda_tangent_2))), (beta_2-(sq((@6*lambda_normal_2))-sq(fabs(lambda_tangent_2)))), (v_0-(p_vt_2-n_vt_2)), (-(v_2+((@4*gamma_3)*lambda_tangent_3))), (beta_3-(sq((@6*lambda_normal_3))-sq(fabs(lambda_tangent_3)))), (v_2-(p_vt_3-n_vt_3)), (-(v_4+((@4*gamma_4)*lambda_tangent_4))), (beta_4-(sq((@6*lambda_normal_4))-sq(fabs(lambda_tangent_4)))), (v_4-(p_vt_4-n_vt_4))]');
            tc.verifyEqual(str(dcs.g_impulse), '@1=-0.707107, @2=0.707107, @3=1, @4=0.5, @5=2, @6=0.1, @7=(sq((@6*Lambda_normal_0))-sq(fabs(Lambda_tangent_0))), @8=(sq((@6*Lambda_normal_1))-sq(fabs(Lambda_tangent_1))), @9=0.2, @10=(sq((@9*Lambda_normal_2))-sq(fabs(Lambda_tangent_2))), @11=(sq((@9*Lambda_normal_3))-sq(fabs(Lambda_tangent_3))), @12=(sq((@9*Lambda_normal_4))-sq(fabs(Lambda_tangent_4))), [(((v_post_impact_0-v_pre_impact_0)-(@1*Lambda_normal_0))-Lambda_tangent_2), (((v_post_impact_1-v_pre_impact_1)-Lambda_normal_2)-((@1*Lambda_tangent_0)+(@1*Lambda_tangent_1))), (((v_post_impact_2-v_pre_impact_2)-((@2*Lambda_normal_0)+(@1*Lambda_normal_1)))-Lambda_tangent_3), (((v_post_impact_3-v_pre_impact_3)-Lambda_normal_3)-(@2*Lambda_tangent_0)), (((v_post_impact_4-v_pre_impact_4)-(@2*Lambda_normal_1))-Lambda_tangent_4), (((v_post_impact_5-v_pre_impact_5)-Lambda_normal_4)-(@2*Lambda_tangent_1)), (Y_gap_0-(((q_2-q_0)-@3)-@3)), (Y_gap_1-(((q_4-q_2)-@3)-@3)), (Y_gap_2-(q_1-@3)), (Y_gap_3-(q_3-@3)), (Y_gap_4-(q_5-@3)), ((P_vn_0-N_vn_0)-((@1*v_post_impact_0)+(@2*v_post_impact_2))), ((P_vn_1-N_vn_1)-((@1*(v_post_impact_2+(@4*v_pre_impact_2)))+(@2*(v_post_impact_4+(@4*v_pre_impact_4))))), ((P_vn_2-N_vn_2)-v_post_impact_1), ((P_vn_3-N_vn_3)-v_post_impact_3), ((P_vn_4-N_vn_4)-v_post_impact_5), (((@2*v_post_impact_1)+(@1*v_post_impact_3))-((@5*Gamma_0)*Lambda_tangent_0)), (Beta_0-@7), (Beta_1-@7), (Beta_2-@7), (Beta_3-@7), (Beta_4-@7), (((@1*v_post_impact_1)+(@2*v_post_impact_3))-(P_vt_0-N_vt_0)), (((@2*v_post_impact_1)+(@1*v_post_impact_5))-((@5*Gamma_1)*Lambda_tangent_1)), (Beta_0-@8), (Beta_1-@8), (Beta_2-@8), (Beta_3-@8), (Beta_4-@8), (((@1*v_post_impact_1)+(@2*v_post_impact_5))-(P_vt_1-N_vt_1)), (-(v_post_impact_0+((@5*Gamma_2)*Lambda_tangent_2))), (Beta_0-@10), (Beta_1-@10), (Beta_2-@10), (Beta_3-@10), (Beta_4-@10), (v_post_impact_0-(P_vt_2-N_vt_2)), (-(v_post_impact_2+((@5*Gamma_3)*Lambda_tangent_3))), (Beta_0-@11), (Beta_1-@11), (Beta_2-@11), (Beta_3-@11), (Beta_4-@11), (v_post_impact_2-(P_vt_3-N_vt_3)), (-(v_post_impact_4+((@5*Gamma_4)*Lambda_tangent_4))), (Beta_0-@12), (Beta_1-@12), (Beta_2-@12), (Beta_3-@12), (Beta_4-@12), (v_post_impact_4-(P_vt_4-N_vt_4))]');
        end
    end

    methods
        function model = simplest_model(tc)
            import casadi.*
            model = nosnoc.model.Cls();

            g = 9.81;
            x0 = [0.8;0];

            q = SX.sym('q',1);
            v = SX.sym('v',1);
            model.M = 1;
            model.x = [q;v];
            model.e = 1;
            model.mu = 0;
            model.x0 = x0;
            model.f_v = -g;
            model.f_c = q;
        end

        function model = carts(tc)
            import casadi.*
            
            m1 = 1;
            m2 = 1;
            m3 = 1;
            cart_width1 = 2;
            cart_width2 = 2;
            cart_width3 = 2;

            M = diag([m1, m1, m2, m2, m3, m3]);

            % Bounds on states and controls
            ubx = ones(12,1)*10;
            lbx = -ones(12,1)*10;
            ubu = 30;
            lbu = -30;

            x0 = [-3; 1; 0; 1;  3; 1; ...
                0; 0; 0; 0; 0; 0];
            u_ref = 0;

            x_ref = [-7; 1; 0; 1; 5; 1;...
                0; 0; 0; 0; 0; 0];
            
            Q = diag([10; 0.001; 1; 0.001; 10; 0.001; ...
                       0.1; 0.1; 0.1; 0.1; 0.1; 0.1]);
            Q_terminal = 100*Q;
            R = 0.1;

            %% Symbolic variables and bounds
            g = 9.81;
            q = SX.sym('q',6);
            v = SX.sym('v',6);
            u = SX.sym('u',1);
            x = [q;v];

            q1 = q(1:2);
            q2 = q(3:4);
            q3 = q(5:6);

            model = nosnoc.model.Cls();
            model.x = x;
            model.u = u;
            model.x0 = x0;

            model.M = M;
            model.f_v = [ 0;...
                -m1*g;
                u;...
                -m2*g;
                0;...
                -m3*g];

            % gap functions
            f_c = [q2(1) - q1(1) - 0.5*cart_width2 - 0.5*cart_width1;...
                q3(1) - q2(1) - 0.5*cart_width3 - 0.5*cart_width2;...
                q1(2)-cart_width1/2;...
                q2(2)-cart_width2/2;...
                q3(2)-cart_width3/2;...
                  ];

            J_tangent = [0  0 1 0 0 ;...
                -1 -1 0 0 0;...
                0  0 0 1 0;...
                1  0 0 0 0;...
                0  0 0 0 1;...
                0  1 0 0 0];

            J_tangent =   J_tangent./vecnorm(J_tangent);
            J_normal = full(f_c.jacobian(q));
            J_normal_fun = Function('J_normal_fun',{q},{J_normal});
            J_normal = full(J_normal_fun(x0(1:6)))';
            J_normal  = J_normal ./vecnorm(J_normal);

            D_tangent  = [];
            for ii = 1:size(J_tangent,2)
                D_tangent  = [D_tangent, J_tangent(:,ii), -J_tangent(:,ii)];
            end

            model.f_c = f_c;
            model.J_normal = J_normal;
            model.J_tangent = J_tangent;
            model.D_tangent = D_tangent;
            model.e =  [0.0 0.5 0.0 0.0 0.0];
            model.mu = [0.1 0.1 0.2 0.2 0.2];
            % box constraints on controls and states
            model.lbu = lbu;
            model.ubu = ubu;
            model.lbx = lbx;
            model.ubx = ubx;
            % Stage cost
            model.f_q = (x-x_ref)'*Q*(x-x_ref)+ u'*R*u;
            model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);

        end
    end
end
