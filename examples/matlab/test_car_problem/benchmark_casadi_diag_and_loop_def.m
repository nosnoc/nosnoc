%% 
import casadi.*

N = 2e3;
x = MX.sym('X',N);

f_1 = diag(x)*x;

f_2 = [];
for ii = 1:N
    f_2 = [f_2;x(ii)*x(ii)];
end

f_1 = Function('f_1',{x},{f_1});
f_2 = Function('f_2',{x},{f_2});


clc


times_1 = [];
times_2 = [];
x0 = rand(N,1);

for jj = 1 :100
tic
f_1_eval = f_1(x0);
times_1 = [times_1;toc];
tic
f_2_eval = f_2(x0);
times_2 = [times_2;toc];
end

close all
figure
stairs(times_1)
grid on
hold on
stairs(times_2)

fprintf('Average 1st: %2.4e 2nd: %2.4e \n ',mean(times_1),mean(times_2))
