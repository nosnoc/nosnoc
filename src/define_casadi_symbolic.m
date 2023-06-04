function sym = define_casadi_symbolic(type, name, dims, sparsity)
    import casadi.*
    if nargin < 3
        dims = 1;
    end

    if ~exist('sparsity', 'var')
        if size(dims, 2) == 1 
            sparsity = Sparsity.dense([dims,1]);
        else
            sparsity = Sparsity.dense(dims);
        end
    end
    
    if strcmp(type, 'casadi.SX') || strcmp(type, 'SX')
        sym = SX.sym(name, sparsity);
    elseif strcmp(type, 'casadi.MX')  || strcmp(type, 'MX')
        sym = MX.sym(name, sparsity);
    else
        error('Type must be MX or SX.')
    end
end

% tests
% type = 'SX'

% x = get_ca_symbolic(type, 'x', [3, 1]);
% disp(x)

% y = get_ca_symbolic(type, 'y', [2, 2]);

% disp(y)

% y = get_ca_symbolic(type, 'y');
% disp(y)
