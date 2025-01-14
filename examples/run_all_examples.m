clear all;
set(0,'DefaultFigureVisible','off');
folder  = '.';
ftext = readlines("all_examples.txt");
[path,name,ext] = fileparts(ftext);
orig_dir = pwd;
%c = parcluster;
for ii=1:length(name)
    cd(path(ii));
    save('ws.mat')
    try
        run(name(ii))
    catch
        load 'ws.mat'
        warning([char(path(ii)) '/' char(name(ii)), ' failed'])
    end
    clear;
    close all;
    %batch(name, 'CaptureDiary', true);
    %c.Jobs
    load 'ws.mat'
    delete 'ws.mat'
    cd(orig_dir);
end

