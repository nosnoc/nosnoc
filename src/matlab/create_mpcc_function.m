% Copyright 2022 Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

% This file is part of NOSNOC.

% The 2-Clause BSD License

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.

% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
function [Psi_mpcc_fun] = create_mpcc_function(mpcc_function,casadi_symbolic_mode)
% Create a CasADi function for the treatment of the bilinear terms in the
% complementarity conditions. It can be the plain bilinear form, some of
% the well known NCP functions, e.g., Fischer-Burmeister and some
% specalized versions, e.g., the Steffenson-Ulbrich local regularization.
%% CasADi
import casadi.*
% %  'FischerBurmeister'; 'NaturalResidual', 'ChenChenKanzow', 'SteffensonUlbrich' , 'KanzowSchwartz' , 'LinFukushima'
%% Variables
a = define_casadi_symbolic(casadi_symbolic_mode,'a',1);
b = define_casadi_symbolic(casadi_symbolic_mode,'b',1);
sigma = define_casadi_symbolic(casadi_symbolic_mode,'sigma',1);

switch mpcc_function
    case 'bilinear'
        Psi_mpcc = a*b-sigma;
    case 'Fischer-Burmeister' 
%         || 'FischerBurmeister'
        Psi_mpcc = a+b-sqrt(a^2+b^2+sigma^2);
%         Psi_mpcc = Psi_mpcc^4;
    case 'Natural-Residual' 
%         || 'min'
        Psi_mpcc = 0.5*(a+b-sqrt((a-b)^2+sigma^2));
%         Psi_mpcc = Psi_mpcc^2;
    case 'Chen-Chen-Kanzow' 
%         || 'ChenChenKanzow'
        alpha = 0.5;
        Psi_mpcc = alpha*(a+b-sqrt(a^2+b^2+sigma^2))+(1-alpha)*(a*b-sigma);

    case 'Steffenson-Ulbrich'
%         || 'SteffensonUlbrich'  || 'Steffenson-Ulbrich-Sin' || 'SteffensonUlbrichSin'
        x = a-b;
        z = x/sigma;
        y_sin = sigma*((2/pi)*sin(z*pi/2+3*pi/2)+1);
        Psi_mpcc = a+b-if_else(abs(x)>=sigma,abs(x),y_sin);
    case  'Steffenson-Ulbrich-Pol'
%         || 'SteffensonUlbrichPol' 
        x = a-b;
        z = x/sigma;
        y_pol = sigma*(1/8*(-z^4+6*z^2+3));
        Psi_mpcc  =a+b- if_else(abs(x)>=sigma,abs(x),y_pol);
    case 'Kanzow-Schwartz'
%         || 'KanzowSchwartz'
        a1 = a-sigma;
        b1 = b-sigma;
        Psi_mpcc  = if_else((a1+b1)>=0,a1*b1,-0.5*(a1^2+b1^2)); 
    case 'Lin-Fukushima' 
%         Psi_mpcc = [a*b-sigma;...
%                     -((a-sigma)*(b-sigma)-sigma^2)];
        Psi_mpcc1 = [a*b-sigma];
        Psi_mpcc2 = [-((a-sigma)*(b-sigma)-sigma^2)];
    case 'Kadrani'
        Psi_mpcc = (a-sigma)*(b-sigma);
    otherwise
        error('Provdied mpcc_function input is not valid.')
end

%% Create CasADi function
try
    Psi_mpcc_fun = Function('Psi_mpcc_fun',{a,b,sigma},{Psi_mpcc});
catch
    Psi_mpcc_fun = Function('Psi_mpcc_fun',{a,b,sigma},{Psi_mpcc1,Psi_mpcc2});
end
end