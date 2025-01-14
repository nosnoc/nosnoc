clear all; close all;
%set(0,'DefaultFigureVisible','off');
folder  = '.';
ftext = readlines("all_examples.txt");
[path,name,ext] = fileparts(ftext);
orig_dir = pwd;
c = parcluster;
for ii=1:length(name)
    cd(path(ii));
    jobs{ii} = batch(name(ii), 'CaptureDiary', true);
    name(ii)
    msg = char(formattedDisplayText(jobs));
    update_msg(msg);

    cd(orig_dir);
end
jobs = [jobs{:}];
while true
    msg = char(formattedDisplayText(jobs));
    update_msg(msg);
    all_done = true;
    for job=jobs
        all_done = all_done & strcmp(job.State, 'finished');
    end
    if all_done
        break
    end
    pause(1);
end

function update_msg(msg)
    ASCII_BKSP_CHAR = 8;
    persistent prev_len;
    if isempty(prev_len) 
        prev_len = 0;
    end
    
    disp([ char(repmat(ASCII_BKSP_CHAR,1,prev_len)) msg]);
    prev_len = numel(msg)+1;
end
