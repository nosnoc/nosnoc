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

function plot_two_gene(results, arrows)
% Warning plotting arrow is _very_ slow because matlab is weird.
    hold on
    for result=results
        plot(result.x_res(1,:), result.x_res(2,:));
        if arrows
            dx = gradient(result.x_res(1,:));
            dy = gradient(result.x_res(2,:));
            quiv = quiver(result.x_res(1,:),result.x_res(2,:),dx,dy,10^-10);

            % Hack to get arrows
            U = quiv.UData;
            V = quiv.VData;
            X = quiv.XData;
            Y = quiv.YData;

            %right version (with annotation)
            head_width = 2.5;
            head_length = 2.5;
            line_length = 1e-10;
            for ii = 1:size(X,1)
                for ij = 1:size(X,2)
                    ah = annotation('arrow',...
                                    'headStyle','cback1','HeadLength',head_length,'HeadWidth',head_width,...
                                    'Color', quiv.Color);
                    set(ah,'parent',gca);
                    set(ah,'position',[X(ii,ij) Y(ii,ij) line_length*U(ii,ij) line_length*V(ii,ij)]);

                end
            end
        end
    end

    xline([4,8],'--',{'$\theta_1^1$','$\theta_1^2$'},'Interpreter','latex')
    yline([4,8],'--',{'$\theta_2^1$','$\theta_2^2$'},'Interpreter','latex')
    xlabel('$x_1$','Interpreter','latex');
    ylabel('$x_2$','Interpreter','latex');
end
