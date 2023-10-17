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

function [results]= form_structured_output(problem, w_opt, name, results)
    ind = problem.(strcat('ind_', name));

    % generate structured output
    opt_s = cellfun(@(idx) w_opt(idx), ind(:,:,:), 'uni', 0);
    opt_s_end = cellfun(@(idx) w_opt(idx), ind(:,:,end), 'uni', 0);
    
    lens = cellfun(@(c) length(c), opt_s);
    if all(lens==0)
        % TODO this is needed for gauss-legendre but causes some problems elsewhere
        %opt_s = cellfun(@(idx) w_opt(idx), ind(:,:,end-1), 'uni', 0);
        return
    end
    opt_s_concat = opt_s(:,:,end);
    for ii=2:size(opt_s,3)
        for jj=1:size(opt_s,1)
            for kk=1:size(opt_s,2)
                opt_s_concat{jj,kk} = horzcat(opt_s_concat{jj,kk}, opt_s{jj,kk,ii});
            end
        end
    end
    results.structured.(name) = opt_s_concat;

    temp = opt_s_end';
    flat = horzcat(temp{:});
    results.(name) = flat;

    % generate extended
    opt_extended = cellfun(@(idx) w_opt(idx), sort_ind_sets(ind(:)), 'uni', 0);

    temp = opt_extended';
    results.extended.(name) = horzcat(temp{:});
end
