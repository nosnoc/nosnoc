function [B,C,D,tau_root] = collocation_times_fesd(n_s,irk_scheme)


import casadi.*

% Degree of interpolating polynomial
% Get collocation points
if isequal(irk_scheme,'lobbato')
    switch n_s
        case 1
            tau_root = [0  1];
        case 2
            tau_root = [0 1/2 1];
        case 3
            tau_root = [0 1/2-sqrt(5)/10 1/2+sqrt(5)/10 1];
        case 4
            tau_root = [0 1/2-sqrt(21)/14  1/2 1/2+sqrt(21)/14 1];
        otherwise
            error('Not implemented')
    end
    irk_scheme = 'radau';
else
tau_root = [0 collocation_points(n_s, irk_scheme)];
end

% Coefficients of the collocation equation
C = zeros(n_s+1,n_s+1);

% Coefficients of the continuity equation
D = zeros(n_s+1, 1);

% Coefficients of the quadrature function
B = zeros(n_s+1, 1);

% Construct polynomial basis
for j=1:n_s+1
    % Construct Lagrange polynomials to get the polynomial basis at the collocation point
    coeff = 1;
    for r=1:n_s+1
        if r ~= j
            coeff = conv(coeff, [1, -tau_root(r)]);
            coeff = coeff / (tau_root(j)-tau_root(r));
        end
    end
    % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D(j) = polyval(coeff, 1.0);
    
    % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    pder = polyder(coeff);
    for r=1:n_s+1
        C(j,r) = polyval(pder, tau_root(r));
    end
    
    % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
    pint = polyint(coeff);
    B(j) = polyval(pint, 1.0);
end

end

