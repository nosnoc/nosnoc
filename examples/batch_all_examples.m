clear all; close all;
%set(0,'DefaultFigureVisible','off');
folder  = '.';
ftext = readlines("all_examples.txt");
[fdir,name,ext] = fileparts(ftext);
orig_dir = pwd;
c = parcluster;
c.NumWorkers = feature('numcores');
for ii=1:2%length(name)
    cd(fdir(ii));
    job = batch(name(ii), 'CaptureDiary', true, 'AutoAttachFiles', false);
    jobs(ii) = job;
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

n_failed = 0;
% Log failures and sucesses
for ii=1:2%length(name)
    if isempty(jobs(ii).Tasks.Error)
        disp([char(ftext(ii)) ' ran without errors']);
    else
        disp([char(ftext(ii)) ' ran with error: ' jobs(ii).Tasks.ErrorIdentifier ' -> ' jobs(ii).Tasks.ErrorMessage]);
        n_failed = n_failed + 1;
    end
end

% Build markdown table for
mkdir ../test-results
md_fid = fopen("../test-results/examples.md", 'w');

if md_fid > 0
    if n_failed
        fprintf(md_fid, '| Example | Success |\n');
        fprintf(md_fid, '| --- | --- |\n');

        for ii=1:length(name)
            if ~isempty(jobs(ii).Tasks.Error)
                fprintf(md_fid, '| %s | %s |\n', ftext(ii), [jobs(ii).Tasks.ErrorIdentifier ' -> ' jobs(ii).Tasks.ErrorMessage])
            end
        end
    else
        fprintf(md_fid, 'All examples ran without error');
    end
    fclose(md_fid);
else
    disp('failed to open file')
end


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
