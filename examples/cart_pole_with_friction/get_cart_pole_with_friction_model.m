% BSD 2-Clause License

% Copyright (c) 2022, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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

function model = get_cart_pole_with_friction_model()
    import casadi.*
    model = NosnocModel();
    model.T = 4;    % Time horizon

    % fixed values
    m1 = 1; % cart
    m2 = 0.1; % link
    g = 9.81;
    link_length = 1;

    % symbolics
    q = SX.sym('q', 2);
    v = SX.sym('v', 2);
    x = [q;v];
    u = SX.sym('u', 1); % control

    % switching function c: cart velocity
    c = v(1);
    % Sign matrix % f_1 for c>0, f_2 for c<0
    S = [1; -1];

    % Inertia matrix
    M = [m1 + m2, m2*link_length*cos(q(2));...
        m2 *link_length*cos(q(2)),  m2*link_length^2];
    % Coriolis force
    C = [0, -m2 * link_length*v(2)*sin(q(2));...
        0,   0];

    % Fixed friction force
    F_friction = 2;

    % all forces = Gravity + Control + Coriolis (+Friction)
    f_all = [0; -m2*g*link_length*sin(x(2))] + [u; 0] - C*v;

    % Dynamics forv > 0
    f_1 = [v;...
            inv(M)*(f_all-[F_friction;0])];
    % Dynamics for v<0
    f_2 = [v;...
            inv(M)*(f_all+[F_friction;0])];
    F = [f_1, f_2];

    % specify initial and end state, cost ref and weight matrix
    x0 = [1; 0/180*pi; 0; 0]; % start downwards
    x_ref = [0; 180/180*pi; 0; 0]; % end upwards
    
    % bounds
    ubx = [5; 240/180*pi; 20; 20];
    lbx = [-0.0; -240/180*pi; -20; -20];
    u_max = 30;
    u_ref = 0;
    
    model.F = F;
    model.c = c;
    model.lbx = lbx;
    model.ubx = ubx;
    model.x = x;
    model.x0 =  x0;
    model.u = u;
    model.lbu = -u_max;
    model.ubu = u_max;
    model.S = S;
    
    % Stage cost
    Q = diag([1; 100; 1; 1]);
    Q_terminal = diag([100; 100; 10; 10]);
    R = 1;
    model.f_q = (x-x_ref)'*Q*(x-x_ref)+ (u-u_ref)'*R*(u-u_ref);
    model.f_q_T = (x-x_ref)'*Q_terminal*(x-x_ref);
end