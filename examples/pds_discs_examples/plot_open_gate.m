function plot_open_gate(h,pos,r,type,fig,vidname)
    n = length(r);
    if ~exist('type')
        type = [];
        for ii=1:n
            type = [type, "circle"];
        end
    end
    if length(type) ~= n
        error('not all types specified')
    end
    indices = {};
    curr_idx = 1;
    for ii=1:n
        switch type(ii)
          case "circle"
            indices{ii} = curr_idx:curr_idx+1;
            curr_idx = curr_idx+2;
          case "square"
            indices{ii} = curr_idx:curr_idx+2;
            curr_idx = curr_idx+3;
        end
    end
    if ~exist('fig')
        fig = figure;
    else
        figure(fig)
    end
    axis equal
    box on;
    ax = gca;
    set(ax,'XTick',[-10,0,10], 'YTick', [0,5,10], 'TickLength', [0.025,0.025], 'LineWidth', 4.0);
    xlim([-12,12])
    ylim([-6,5])
    ax.XAxis.FontSize = 36;
    ax.YAxis.FontSize = 36;
    xlabel("$c_x$", "fontsize", 32)
    ylabel("$c_y$", "fontsize", 32)
    
    %axis off
    discs = {};
    tau = linspace(0, 2*pi)';
    rot = @(theta) [cos(theta) -sin(theta);...
        sin(theta) cos(theta)];
    circ = @(p,r) [r*cos(tau)+p(1),r*sin(tau)+p(2)];
    sqr = @(p,r) transpose(rot(-p(3))*([r,r;r,-r;-r,-r;-r,r]') + p(1:2));
    for ii=1:n
        p = pos(indices{ii},1);
        switch type(ii)
          case "circle"
            v = circ(p,r(ii));
          case "square"
            v = sqr(p,r(ii));
        end
        if ii==3
            color = [0.8500 0.3250 0.0980];
        else
            color = [0 0.4470 0.7410];
        end
        discs{ii} = patch(v(:,1),v(:,2), color);
    end
    % plot "gate"
    gate = patch([-20,20,20,-20],[-10,-10,pos(end,1), pos(end,1)], [1 0 0], 'FaceAlpha', 0.1, 'LineStyle', '--', 'EdgeColor' , [1 0 0]);
    
    
    if exist('vidname')
        mkdir([vidname '_frames'])
        writer = VideoWriter([vidname '.avi']);
        open(writer);
        frame = getframe(gca);
        exportgraphics(gca, [vidname '_frames/' num2str(1) '.pdf'])
        writeVideo(writer,frame)
    end
    pause(h(1));
    for jj=2:(length(h))
        for ii=1:n
            p = pos(indices{ii},jj+1);
            switch type(ii)
              case "circle"
                v = circ(p,r(ii));
              case "square"
                v = sqr(p,r(ii));
            end
            discs{ii}.Vertices = v;
        end
        gate.Vertices = [-20 -10;
            20 -10;
            20 pos(end,jj);
            -20 pos(end,jj)];
        drawnow limitrate;
        %pause(h(jj));
        pause(0.1)
        if exist('vidname')
            frame = getframe(gca);
            ff=jj+1;
            %exportgraphics(gca, [vidname '_frames/' num2str(jj) '.pdf'])
            writeVideo(writer,frame);
        end
    end
    if exist('vidname')
        close(writer);
    end
end
