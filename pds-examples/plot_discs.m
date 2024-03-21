function plot_discs(h,pos,r,type,lim,fig,vidname)
    n = length(r);
    if ~exist('lim')
        lim = 50;
    end
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
    set(gca,'XTick',[], 'YTick', [])
    xlim([-lim,lim])
    ylim([-lim,lim])
    axis off
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
        discs{ii} = patch(v(:,1),v(:,2),'g');
    end
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
        drawnow limitrate;
        %pause(h(jj));
        pause(0.1)
        if exist('vidname')
            frame = getframe(gca);
            ff=jj+1;
            exportgraphics(gca, [vidname '_frames/' num2str(jj) '.pdf'])
            writeVideo(writer,frame);
        end
    end
    if exist('vidname')
        close(writer);
    end
end
