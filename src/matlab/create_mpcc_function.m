%    This file is part of NOSNOC.
%
%    NOSNOC -- A software for NOnSmooth Numerical Optimal Control.
%    Copyright (C) 2022 Armin Nurkanovic, Moritz Diehl (ALU Freiburg).
%
%    NOSNOC is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 3 of the License, or (at your option) any later version.
%
%    NO-NOC is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with NOS-NOC; if not, write to the Free Software Foundation,
%    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
function [Psi_mpcc_fun] = create_mpcc_function(mpcc_function,casadi_symbolic_mode,sigma_p)
% Create a CasADi function for the treatment of the bilinear terms in the
% complementarity conditions. It can be the plain bilinaer form, some of
% the well known NCP functions, e.g., Fischer-Burmeister and some
% specalized versions, e.g., the Steffenson-Ulbrich local regularization.
%% CasADi
import casadi.*
% %  'FischerBurmeister'; 'NaturalResidual', 'ChenChenKanzow', 'SteffensonUlbrich' , 'KanzowSchwartz' , 'LinFukushima'
%% Variables
a = define_casadi_symbolic(casadi_symbolic_mode,'a',1);
b = define_casadi_symbolic(casadi_symbolic_mode,'b',1);

switch mpcc_function
    case 'bilinear'
        Psi_mpcc = a*b-sigma_p;
    case 'Fischer-Burmeister' 
%         || 'FischerBurmeister'
        Psi_mpcc = a+b-sqrt(a^2+b^2-sigma_p^2);
    case 'Natural-Residual' 
%         || 'min'
        Psi_mpcc = 0.5*(a+b-sqrt((a-b)^2+sigma_p^2));
    case 'Chen-Chen-Kanzow' 
%         || 'ChenChenKanzow'
        alpha = 0.5;
        Psi_mpcc = alpha*(a+b-sqrt(a^2+b^2-sigma_p^2))+(1-alpha)*(a*b-sigma_p);
    case 'Steffenson-Ulbrich'
%         || 'SteffensonUlbrich'  || 'Steffenson-Ulbrich-Sin' || 'SteffensonUlbrichSin'
        x = a-b;
        z = x/sigma_p;
        y_sin = sigma_p*((2/pi)*sin(z*pi/2+3*pi/2)+1);
        Psi_mpcc = if_else(abs(x)>=sigma,abs(x),y_sin);
    case  'Steffenson-Ulbrich-Pol'
%         || 'SteffensonUlbrichPol' 
        x = a-b;
        z = x/sigma_p;
        y_pol = sigma_p*(1/8*(-z^4+6*z^2+3));
        Psi_mpcc  = if_else(abs(x)>=sigma,abs(x),y_pol);
    case 'Kanzow-Schwartz'
%         || 'KanzowSchwartz'
        %
    case 'Lin-Fukushima' 
%         || 'LinFukushima'
        %
    case 'Kadrani'
        %
    otherwise
        error('Provdied mpcc_function input is not valid.')
end

%% Create CasADi function
Psi_mpcc_fun = Function('Psi_mpcc_fun',{a,b,sigma_p},{Psi_mpcc});
end