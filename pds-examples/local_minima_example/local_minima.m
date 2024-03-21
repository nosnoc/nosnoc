%clear all
close all
import vdx.*
N_sim = 1;
ts = 0;
xs = 1.0:0.01:2.0;%-1.6:0.1:-1;
x_term = [];

f1 = figure;
legend('location', 'north')
f2 = figure;
%fplot(@(x) 1/24*(256 + 96*x + 12*x.^2 + x.^3) - (x.^3)/8 + (4+x/2 - 3.5).^2, "DisplayName", 'analytical')
%fplot(@(x) (4+x/2 - 3.5).^2, "DisplayName", 'analytical')
%xlim([-1.5 -0.5])
legend('location', 'north')

for comp_tol=[1e-6]%,1e-3,1e-6,1e-9]
    xplot = [1];
    fs = [1];
    ress = [1];
    figure(f1)
    p1 = plot(xplot, fs, "DisplayName", ['$\sigma =' num2str(comp_tol) '$']);
    p1.XDataSource = 'xplot';
    p1.YDataSource = 'fs';
    figure(f2)
    p2 = plot(xplot, ress, "DisplayName", ['$\sigma =' num2str(comp_tol) '$']);
    p2.XDataSource = 'xplot';
    p2.YDataSource = 'ress';
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
            xplot(ii) = x0;
            x_term = [x_term,x_curr];
            %f = cost + (x_curr(1)-3.5)^2;
            %disp(x_curr)
            f = (x_curr(2)-1.2)^2;
            disp(x_res)
            ress(ii) = x_curr(2);
            fs(ii) = f;
            figure(f1)
            refreshdata
            figure(f2)
            refreshdata
            %drawnow
            ii=ii+1;
        end
    end
    
end
