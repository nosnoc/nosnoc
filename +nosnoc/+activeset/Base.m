classdef Base < handle
% Base class for all active set objects, currently used just for type verification.

% Propreties of all other activeset classes have different names because they are different. 
% I.e., in PDS/CLS it is active unilateral constraints that define the active set, 
% wheras in PSS it is specifically the active regions, and in Heaviside case it is 
% which of the three branches of the Heaviside step each step function is on.

% TODO(@anton): make this non-instantiable perhaps with a protected constructor
end
