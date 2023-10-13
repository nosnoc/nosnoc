% BSD 2-Clause License

% Copyright (c) 2023, Armin Nurkanović, Jonathan Frey, Anton Pozharskiy, Moritz Diehl

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
classdef CFunctionType
    enumeration
        SCHOLTES
        SCHOLTES_TWO_SIDED
        FISCHER_BURMEISTER
        NATURAL_RESIDUAL
        CHEN_CHEN_KANZOW
        STEFFENSEN_ULBRICH
        STEFFENSEN_ULBRICH_POLY
        KANZOW_SCHWARTZ
        LIN_FUKUSHIMA
        KADRANI
    end

    methods
        function latex = to_latex(obj)
            switch obj
                case CFunctionType.SCHOLTES
                    latex = '$\phi_{\textrm{scholtes}}$';
                case CFunctionType.SCHOLTES_TWO_SIDED
                    latex = '$\phi_{\textrm{two-sided}}$';
                case CFunctionType.FISCHER_BURMEISTER
                    latex = '$\phi_{\textrm{fb}}$';
                case CFunctionType.NATURAL_RESIDUAL
                    latex = '$\phi_{\textrm{nr}}$';
                case CFunctionType.CHEN_CHEN_KANZOW
                    latex = '$\phi_{\textrm{cck}}$';
                case CFunctionType.STEFFENSEN_ULBRICH
                    latex = '$\phi_{\textrm{su1}}$';
                case CFunctionType.STEFFENSEN_ULBRICH_POLY
                    latex = '$\phi_{\textrm{su2}}$';
                case CFunctionType.KANZOW_SCHWARTZ
                    latex = '$\phi_{\textrm{ks}}$';
                case CFunctionType.LIN_FUKUSHIMA
                    latex = '$\phi_{\textrm{lf}}$';
                case CFunctionType.KADRANI
                    latex = '$\phi_{\textrm{ns}}$';
            end
        end

        function acronym = acronym(obj)
            switch obj
                case CFunctionType.SCHOLTES
                    acronym = 'Scholtes';
                case CFunctionType.SCHOLTES_TWO_SIDED
                    acronym = 'Scholtes 2-sided$';
                case CFunctionType.FISCHER_BURMEISTER
                    acronym = 'FB';
                case CFunctionType.NATURAL_RESIDUAL
                    acronym = 'NR';
                case CFunctionType.CHEN_CHEN_KANZOW
                    acronym = 'CCK';
                case CFunctionType.STEFFENSEN_ULBRICH
                    acronym = 'SU1';
                case CFunctionType.STEFFENSEN_ULBRICH_POLY
                    acronym = 'SU2';
                case CFunctionType.KANZOW_SCHWARTZ
                    acronym = 'KS';
                case CFunctionType.LIN_FUKUSHIMA
                    acronym = 'LF';
                case CFunctionType.KADRANI
                    acronym = 'Kadrani';
            end
        end

        function full = full(obj)
            switch obj
                case CFunctionType.SCHOLTES
                    full = 'Scholtes';
                case CFunctionType.SCHOLTES_TWO_SIDED
                    full = 'Scholtes 2-sided';
                case CFunctionType.FISCHER_BURMEISTER
                    full = 'Fischer-Burmeister';
                case CFunctionType.NATURAL_RESIDUAL
                    full = 'Natural Residual';
                case CFunctionType.CHEN_CHEN_KANZOW
                    full = 'Chen-Chen-Kanzow';
                case CFunctionType.STEFFENSEN_ULBRICH
                    full = 'Steffensen-Ulbrich Trigonometric';
                case CFunctionType.STEFFENSEN_ULBRICH_POLY
                    full = 'Steffensen-Ulbrich Polynomial';
                case CFunctionType.KANZOW_SCHWARTZ
                    full = 'Kanzow-Schwartz';
                case CFunctionType.LIN_FUKUSHIMA
                    full = 'Lin-Fukushima';
                case CFunctionType.KADRANI
                    full = 'Kadrani';
            end
        end
    end
end
