%clear all
close all
import vdx.*
N_sim = 10;
ts = 0;
xs = 0.1:0.1:3.5;%-1.6:0.1:-1;
x_term = [];
tols = [1e-3,1e-6,1e-9];
f1 = figure;
xlabel("$x_2(0)$")
ylabel("$x_2(1;x_1(0))$")
title("Nonsmooth Sensitivity")
%ylim([.95,1.3])
xlim([0.1,3.5])

legend('location', 'north')
jj = 1;
xplot = {};
x2_T = {};

for n_s = [2]
    for comp_tol=tols
        hold on
        xplot{jj} = [2.5];
        x2_T{jj} = [1.25];
        p1 = plot(xplot{1}, x2_T{1}, "DisplayName", ['$\sigma =' num2str(comp_tol) '$ $n_s = ' num2str(n_s) '$']);
        p1.XDataSource = ['xplot{' num2str(jj) '}'];
        p1.YDataSource = ['x2_T{' num2str(jj) '}'];
        xplot{jj} = [];
        x2_T{jj} = [];
        ii = 1
        for x0=xs
            for t0=ts
                %disp(x0)
                [prob,data,opts,h] = integrator2(false,N_sim, x0, n_s, false);
                x_curr = [x0;0];
                x_res = x_curr;
                lambda_curr = 0;
                cost = 0;
                for step=1:N_sim
                    prob.w.x(0,0,data.n_s).init = x_curr;
                    prob.w.x(0,0,data.n_s).lb = x_curr;
                    prob.w.x(0,0,data.n_s).ub = x_curr;
                    prob.w.lambda(0,0,data.n_s).init = lambda_curr;
                    prob.w.lambda(0,0,data.n_s).lb = lambda_curr;
                    prob.w.lambda(0,0,data.n_s).ub = lambda_curr;
                    success = homotopy(prob, comp_tol, comp_tol);
                    %success = homotopy(prob);
                    if ~success
                        disp(['Failure to converge at step=' num2str(step)])
                    else
                        %disp(['step=' num2str(step)])
                    end
                    x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
                    lambda_curr = prob.w.lambda(1,data.N_fe,data.n_s).res;
                    x_sim = prob.w.x(1,:,data.n_s).res;
                    x_res = [x_res,x_sim];
                    %disp(prob.g.objective(1).res)
                    %cost = cost + prob.f_result;
                    %cost = cost + prob.g.objective(1).res;
                end
                xplot{jj} = [xplot{jj}, x0];
                x2_T{jj} = [x2_T{jj},x_curr(2)];
                x_term = [x_term,x_curr];
                %f = cost + (x_curr(1)-3.5)^2;
                disp(x_res)
                %disp(prob.w.lambda.res)
                refreshdata
                %drawnow
                ii=ii+1;
            end
        end
        jj = jj+1;
    end
end


hold on
xplot{jj} = [2.5];
x2_T{jj} = [1.25];
p1 = plot(xplot{1}, x2_T{1}, "DisplayName", ['FESD']);
p1.XDataSource = ['xplot{' num2str(jj) '}'];
p1.YDataSource = ['x2_T{' num2str(jj) '}'];
xplot{jj} = [];
x2_T{jj} = [];
ii = 1
for x0=[1]
    for t0=ts
        %disp(x0)
        [prob,data,opts,h] = integrator2(true,N_sim, x0, 2, false);
        x_curr = [x0;0];
        x_res = x_curr;
        h_res = [];
        lambda_curr = 0;
        cost = 0;
        prob.p.gamma_h(1).init = 1e3;
        %prob.w.h(1).init = 0.02;
        %prob.w.h(2).init = 1-0.02;
        for step=1:N_sim
            prob.w.x(0,0,data.n_s).init = x_curr;
            prob.w.x(0,0,data.n_s).lb = x_curr;
            prob.w.x(0,0,data.n_s).ub = x_curr;
            prob.w.lambda(0,0,data.n_s).init = lambda_curr;
            prob.w.lambda(0,0,data.n_s).lb = lambda_curr;
            prob.w.lambda(0,0,data.n_s).ub = lambda_curr;
            success = homotopy(prob,1, 1e-12);
            %success = homotopy(prob);
            if ~success
                disp(['Failure to converge at step=' num2str(step)])
            else
                %disp(['step=' num2str(step)])
            end
            x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
            lambda_curr = prob.w.lambda(1,data.N_fe,data.n_s).res;
            x_sim = prob.w.x(1,:,data.n_s).res;
            x_res = [x_res,x_sim];
            h_res = [h_res,prob.w.h.res];
            %disp(prob.g.objective(1).res)
            %cost = cost + prob.f_result;
            %cost = cost + prob.g.objective(1).res;
        end
        xplot{jj} = [xplot{jj}, x0];
        x2_T{jj} = [x2_T{jj},x_curr(2)];
        x_term = [x_term,x_curr];
        %f = cost + (x_curr(1)-3.5)^2;
        %disp(x_curr)
        disp(prob.w.lambda.res)
        refreshdata
        %drawnow
        ii=ii+1;
    end
end
t_res = [0,cumsum(h_res)];
