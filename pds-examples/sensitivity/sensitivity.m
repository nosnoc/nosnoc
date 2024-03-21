%clear all
close all
import vdx.*
N_sim = 1;
ts = 0;
xs = 2.5:0.01:3.5;%-1.6:0.1:-1;
x_term = [];
tols = [1e-1,5e-2,1e-2,5e-3,1e-3,1e-6,1e-9];
f1 = figure;
xlabel("$x_2(0)$")
ylabel("$x_2(1;x_1(0))$")
title("Nonsmooth Sensitivity")
ylim([.95,1.3])
xlim([2.5,3.5])

legend('location', 'north')
jj = 1;
xplot = {};
x2_T = {};

for comp_tol=tols
    hold on
    xplot{jj} = [2.5];
    x2_T{jj} = [1.25];
    p1 = plot(xplot{1}, x2_T{1}, "DisplayName", ['$\sigma =' num2str(comp_tol) '$']);
    p1.XDataSource = ['xplot{' num2str(jj) '}'];
    p1.YDataSource = ['x2_T{' num2str(jj) '}'];
    xplot{jj} = [];
    x2_T{jj} = [];
    ii = 1
    for x0=xs
        for t0=ts
            %disp(x0)
            [prob,data,opts,h] = integrator(false,N_sim, x0);
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
            %disp(x_curr)
            disp(x_res)
            refreshdata
            %drawnow
            ii=ii+1;
        end
    end
    jj = jj+1;
end

% analytical
analytical = interp1([2.5,3,3.5],[1.25,1,1], xplot{1})

plot(xplot{1}, analytical, 'DisplayName', 'analytical');


figure
legend('location', 'northeastoutside')
jj = 1;
hold on
for comp_tol=tols
    delta = x2_T{jj} - analytical;
    plot(xplot{1}, delta, "DisplayName", ['$\sigma =' num2str(comp_tol) '$'])
    jj = jj+1;
end
