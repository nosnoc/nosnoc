clear all; close all;
%set(0,'DefaultFigureVisible','off');
folder  = '.';
ftext = readlines("all_examples.txt");
[path,name,ext] = fileparts(ftext);
orig_dir = pwd;
c = parcluster;
for ii=1:5%length(name)
    cd(path(ii));
    job = batch(name(ii), 'CaptureDiary', true);
    jobs(ii) = job;
    name(ii)
    %msg = char(formattedDisplayText(jobs, 'SuppressMarkup', true));
    %update_msg(msg);

    cd(orig_dir);
end

while true
    msg = char(formattedDisplayText(jobs, 'SuppressMarkup', true));
    update_msg(msg);
    all_done = true;
    for job=jobs
        all_done = all_done & strcmp(job.State, 'finished');
    end
    if all_done
        break
    end
    pause(10);
end

% Log failures and sucesses
for ii=1:5%length(name)
    if isempty(jobs(ii).Tasks.Error)
        disp([char(ftext(ii)) ' ran without errors']);
    else
        disp([char(ftext(ii)) ' ran with error: ' jobs(ii).Tasks.Error.identifier ' ->' jobs(ii).Tasks.Error.message]);
    end
end

% Build markdown table for
md_fid = fopen("../test-results/examples.md");
jobs(1)

function update_msg(msg)
    ASCII_BKSP_CHAR = 8;
    persistent prev_len;
    if isempty(prev_len) 
        prev_len = 0;
    end
    
    %disp([ char(repmat(ASCII_BKSP_CHAR,1,prev_len)) msg]);
    disp(msg);
    prev_len = numel(msg)+1;
end
