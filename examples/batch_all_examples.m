clear all; close all;
%set(0,'DefaultFigureVisible','off');
folder  = '.';
ftext = readlines("all_examples.txt");
[path,name,ext] = fileparts(ftext);
orig_dir = pwd;
c = parcluster;
for ii=1:length(name)
    cd(path(ii));
    jobs(ii) = batch(name(ii), 'CaptureDiary', true);
    cd(orig_dir);
end

while true
    all_done = true;
    for job=jobs
        all_done = all_done & strcmp(job.State, 'finished');
    end
    if all_done
        break
    end
    pause(1);
    jobs.display
end
