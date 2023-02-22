function plot_circle(varargin)
x = varargin{1};
y = varargin{2};
r = varargin{3};
if nargin > 3
    color = varargin{4};
end


hold on

for ii = 1:length(x)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x(ii);
    yunit = r * sin(th) + y(ii);
    if nargin > 3
        h = plot(xunit, yunit,'color',color);
%         patch([xunit*1e-3 flip(xunit)],[yunit*1e-3 flip(yunit)],'k')
    else
        h = plot(xunit, yunit,'LineWidth',2);
    end
end
end