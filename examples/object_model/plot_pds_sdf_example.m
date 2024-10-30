function plot_pds_sdf_example(h,pos,p_d,indices,pgons,facecolors,linecolors,fig,vidname,export_pdf)
% TODO(@anton) use argument parser
    n_interp = 5;
    m = length(pos);
    t = 1:m;
    t_new = 1:(1/n_interp):m;
    m_new = length(t_new);

    pos = interp1(t, pos', t_new)';
    p_d = interp1(t, p_d', t_new)';
    n = length(pgons);
    if ~exist('export_pdf')
        export_pdf = false;
    end
    if ~exist('fig')
        fig = figure;
    else
        figure(fig)
    end
    axis equal
    box on;
    ax = gca;
    set(ax,'XTick',[0,5,10], 'YTick', [0,5,10], 'TickLength', [0.025,0.025], 'LineWidth', 4.0);
    xlim([-4,4])
    ylim([-4,4])
    % set(ax,'XTick',[-2.5,0,2.5], 'YTick', [-2.5,0,2.5], 'TickLength', [0.025,0.025], 'LineWidth', 4.0);
    % xlim([-5,5])
    % ylim([-5,5])
    ax.XAxis.FontSize = 52;
    ax.YAxis.FontSize = 52;
    xlabel("$x$", "fontsize", 52)
    ylabel("$y$", "fontsize", 52)

    
    hold on
    for ii=1:n
        p = pos(indices{ii},1);
        switch length(indices{ii})
          case 2
            pgonplots{ii} = plot(translate(pgons(ii), p'), 'FaceColor', facecolors{ii}, 'FaceAlpha', 1, 'EdgeColor', linecolors{ii});
          case 3
            pgonplots{ii} = plot(translate(rotate(pgons(ii),rad2deg(p(3))), p(1:2)'), 'FaceColor', facecolors{ii}, 'FaceAlpha', 1, 'EdgeColor', linecolors{ii});
        end
    end
    if ~isempty(p_d)
        pd_line = plot(p_d(1:2:end,1), p_d(2:2:end,1), 'r.', 'markersize', 40);
        pd_line.XDataSource = 'p_d(1:2:end,jj+1)';
        pd_line.YDataSource = 'p_d(2:2:end,jj+1)';
    end
    hold off
    if exist('vidname')
        mkdir([vidname '_frames'])
        writer = VideoWriter([vidname '.avi']);
        writer.Quality = 100;
        writer.FrameRate = 8*n_interp/(mean(h));
        open(writer);
        frame = getframe(gca);
        if export_pdf
            exportgraphics(gca, [vidname '_frames/' num2str(1) '.pdf'])
        end
        writeVideo(writer,frame)
    end
    pause(h(1));
    for jj=1:m_new-1
        for ii=1:n
            p = pos(indices{ii},jj+1);
            switch length(indices{ii})
              case 2
                pgonplots{ii}.Shape = translate(pgons(ii), p');
              case 3
                pgonplots{ii}.Shape = translate(rotate(pgons(ii),rad2deg(p(3))), p(1:2)');
            end
        end
        drawnow limitrate;
        if ~isempty(p_d)
            refreshdata(pd_line,'caller');
        end
        pause(0.001);
        %pause(0.01)
        if exist('vidname')
            frame = getframe(gca);
            ff=jj+1;
            if export_pdf
                set(gca,'Units','normalized','OuterPosition',[0 0 1 1]);
                exportgraphics(gca, [vidname '_frames/' vidname num2str(jj) '.pdf'])
            end
            writeVideo(writer,frame);
        end
    end
    if exist('vidname')
        close(writer);
    end
end
