% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% This file is part of NOSNOC.
classdef NosnocSolverOptions < handle
% TODO clean up much of the work here.
    properties
        % General
        solver_name {mustBeTextScalar} = 'nosnoc_solver'
        nlpsol = 'ipopt'
        casadi_symbolic_mode {mustBeMember(casadi_symbolic_mode,{'casadi.SX', 'casadi.MX'})} = 'casadi.SX'

        % MPCC and Homotopy Settings
        comp_tol(1,1) double {mustBeReal, mustBePositive} = 1e-9
        mpcc_mode(1,1) MpccMode = MpccMode.Scholtes_ineq % 'direct', 'Scholtes_eq', 'Scholtes_ineq', 'ell_1_penalty', 'elastic_ineq', 'elastic_eq' , 'elastic_two_sided',
                                           % 'elastic_ell_1_ineq', 'elastic_ell_1_eq', 'elastic_ell_1_two_sided'
        objective_scaling_direct(1,1) logical = 1
        sigma_0(1,1) double {mustBeReal, mustBePositive} = 1
        sigma_N(1,1) double {mustBeReal, mustBePositive} = 1e-9
        homotopy_update_rule = 'linear' % 'linear' sigma_k = homotopy_update_slope*sigma_N
                                        % 'superlinear' - sigma_k = max(sigma_N,min(homotopy_update_slope*sigma_k,sigma_k^homotopy_update_exponent))
                                        % TODO enum
        homotopy_update_slope(1,1) double {mustBeReal, mustBeInRange(homotopy_update_slope, 0, 1, 'exclusive')} = 0.1
        homotopy_update_exponent(1,1) double {mustBeReal, mustBePositive} = 1.5 % the exponent in the superlinear rule
        N_homotopy = 0 % 0 -> set automatically
        s_elastic_max(1,1) double {mustBeReal, mustBePositive} = 1e1
        s_elastic_min(1,1) double {mustBeReal, mustBeNonnegative} = 0
        s_elastic_0(1,1) double {mustBeReal, mustBePositive} = 1
        elastic_scholtes(1,1) logical = 0

        polishing_step(1,1) logical = 0 % heuristic for fixing active set, yet exerimental, not recommended to use.
        polishing_derivative_test(1,1) logical = 0 % check in sliding mode also the derivative of switching functions

        % Verbose
        print_level = 3
        print_details_if_infeasible = 0;
        pause_homotopy_solver_if_infeasible = 0;

        % IPOPT Settings
        tol_ipopt(1,1) double {mustBeReal, mustBePositive} = 1e-16
        opts_casadi_nlp

        % Integrator Specific
        use_previous_solution_as_initial_guess(1,1) logical = 0
        simulation_problem(1,1) logical = 0

        % TODO: make proper multiple solver class.
        multiple_solvers(1,1) logical = 0

        % All NLP parameters
        p_val

        % psi func
        psi_fun_type CFunctionType = CFunctionType.BILINEAR
        relaxation_method(1,1) RelaxationMode = RelaxationMode.INEQ
        elasticity_mode(1,1) ElasticityMode = ElasticityMode.NONE
        psi_fun

        % Output options
        store_integrator_step_results(1,1) logical = 0
        ipopt_callback = [] % This should be a function handle that takes (model,problem,settings,ipopt_solver,results)
    end

    properties(Dependent)
        time_rescaling
        
    end

    methods
        function obj = NosnocProblemOptions()

            default_tol = 1e-12;

            obj.opts_casadi_nlp.ipopt.print_level = 0;
            obj.opts_casadi_nlp.print_time = 0;
            obj.opts_casadi_nlp.ipopt.sb = 'yes';
            obj.opts_casadi_nlp.verbose = false;
            obj.opts_casadi_nlp.ipopt.max_iter = 500;
            obj.opts_casadi_nlp.ipopt.bound_relax_factor = 0;
            obj.opts_casadi_nlp.ipopt.tol = default_tol;
            obj.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
            obj.opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
            obj.opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
            obj.opts_casadi_nlp.snopt = struct();

            obj.p_val = [obj.sigma_0,obj.rho_sot,obj.rho_h,obj.rho_terminal,obj.T_val];
        end

        function [] = preprocess(obj)
            import casadi.*
            % automatically set up casadi opts
            if obj.print_level < 4 
                obj.opts_casadi_nlp.ipopt.print_level=0;
                obj.opts_casadi_nlp.print_time=0;
                obj.opts_casadi_nlp.ipopt.sb= 'yes';
            elseif print_level == 4
                obj.opts_casadi_nlp.ipopt.print_level=0;
                obj.opts_casadi_nlp.print_time=1;
                obj.opts_casadi_nlp.ipopt.sb= 'no';
            else
                obj.opts_casadi_nlp.ipopt.print_level = 5;
            end

            % MPCC mode setup
            if obj.mpcc_mode == MpccMode.direct
                % TODO maybe need to check later if people try to set sigmas/n_homotopy as for now we just let them break the
                %      assumptions.
                if obj.print_level >= 1
                    fprintf('Info: Setting N_homotopy to 1 and sigma_0,sigma_N to constants in direct mode.\n')
                end
                obj.N_homotopy = 1;
                obj.sigma_0 = 1e-12;
                obj.sigma_N = 1e-12;
                obj.mpcc_mode = MpccMode.Scholtes_ineq;
            end
            if obj.mpcc_mode == MpccMode.ell_1_penalty
                if obj.print_level >= 1
                    fprintf('Info: Setting cross_comp_mode to 12 in ell_1_penalty mode.\n')
                end
                obj.cross_comp_mode = 12;
            end

            switch obj.mpcc_mode
                case MpccMode.Scholtes_ineq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.INEQ;
                    obj.elasticity_mode = ElasticityMode.NONE;
                case MpccMode.Scholtes_eq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.EQ;
                    obj.elasticity_mode = ElasticityMode.NONE;
                case MpccMode.elastic_ineq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.INEQ;
                    obj.elasticity_mode = ElasticityMode.ELL_INF;
                case MpccMode.elastic_eq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.EQ;
                    obj.elasticity_mode = ElasticityMode.ELL_INF;
                case MpccMode.elastic_two_sided
                    obj.psi_fun_type = CFunctionType.BILINEAR_TWO_SIDED;
                    obj.relaxation_method = RelaxationMode.TWO_SIDED;
                    obj.elasticity_mode = ElasticityMode.ELL_INF;
                case MpccMode.elastic_ell_1_ineq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.INEQ;
                    obj.elasticity_mode = ElasticityMode.ELL_1;
                case MpccMode.elastic_ell_1_eq
                    obj.psi_fun_type = CFunctionType.BILINEAR;
                    obj.relaxation_method = RelaxationMode.EQ;
                    obj.elasticity_mode = ElasticityMode.ELL_1;
                case MpccMode.elastic_ell_1_two_sided
                    obj.psi_fun_type = CFunctionType.BILINEAR_TWO_SIDED;
                    obj.relaxation_method = RelaxationMode.TWO_SIDED;
                    obj.elasticity_mode = ElasticityMode.ELL_1;
            end

            % psi function
            a = define_casadi_symbolic(obj.casadi_symbolic_mode,'a',1);
            b = define_casadi_symbolic(obj.casadi_symbolic_mode,'b',1);
            sigma = define_casadi_symbolic(obj.casadi_symbolic_mode,'sigma',1);

            switch obj.psi_fun_type
              case CFunctionType.BILINEAR
                psi_mpcc = a.*b-sigma;
                
              case CFunctionType.BILINEAR_TWO_SIDED
                psi_mpcc = [a*b-sigma,a*b+sigma];
                
              case CFunctionType.FISCHER_BURMEISTER
                psi_mpcc = a+b-sqrt(a^2+b^2+sigma^2);

              case CFunctionType.NATURAL_RESIDUAL
                psi_mpcc = 0.5*(a+b-sqrt((a-b)^2+sigma^2));

              case CFunctionType.CHEN_CHEN_KANZOW
                alpha = 0.5;
                psi_mpcc = alpha*(a+b-sqrt(a^2+b^2+sigma^2))+(1-alpha)*(a*b-sigma);

              case CFunctionType.STEFFENSON_ULBRICH
                x = a-b;
                z = x/sigma;
                y_sin = sigma*((2/pi)*sin(z*pi/2+3*pi/2)+1);
                psi_mpcc = a+b-if_else(abs(x)>=sigma,abs(x),y_sin);

              case  CFunctionType.STEFFENSON_ULBRICH_POLY
                x = a-b;
                z = x/sigma;
                y_pol = sigma*(1/8*(-z^4+6*z^2+3));
                psi_mpcc  =a+b- if_else(abs(x)>=sigma,abs(x),y_pol);

              case CFunctionType.KANZOW_SCHWARTZ
                a1 = a-sigma;
                b1 = b-sigma;
                psi_mpcc  = if_else((a1+b1)>=0,a1*b1,-0.5*(a1^2+b1^2));

              case CFunctionType.LIN_FUKUSHIMAgs
                psi_mpcc1 = [a*b-sigma];
                psi_mpcc2 = [-((a-sigma)*(b-sigma)-sigma^2)];
                psi_mpcc = vertcat(psi_mpcc1, psi_mpcc2)

              case CFunctionType.KADRANI
                psi_mpcc = (a-sigma)*(b-sigma);
            end

            obj.psi_fun = Function('psi_fun',{a,b,sigma},{psi_mpcc});
        end

    end
end