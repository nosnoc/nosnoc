classdef IntegratorType
% Enumerated type for the types of available integrators.
%
% Note:
%       Not every option is valid for every model type. For example `SMOOTHED_PSS` only works for piecewise smooth systems.
    enumeration
        FESD          % Use an FESD discretization and an MPEC/MPCC solver for each integration step. 
        SMOOTHED_PSS  % Smooth a piecewise smooth system with heaviside step functions.
        STEWART       % Stewart's event-based integration method for Filippov PSS systems. 
    end
end
