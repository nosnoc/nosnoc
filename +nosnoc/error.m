function error(error_id, msg)
% nosnoc specific wrapper for error which builds the full error_id automagically.
% This may break if for some reason it is called outside of `+nosnoc` or its subdirectories
% or if the user has a +nosnoc directory on the path somehow (Though this is unlikely due
% to how we walk down the path).

    error_struct.message = ['nosnoc: ' char(msg)];
    
    stack = dbstack('-completenames', 1);

    [last_frame_path,last_frame_fname,~] = fileparts(stack(1).file);

    last_frame_path_parts = split(last_frame_path, filesep);

    id_complete = false;
    error_id = [normalize_part_for_id(last_frame_fname) ':' normalize_part_for_id(error_id)];
    ii = length(last_frame_path_parts);
    while ~id_complete
        if strcmp(last_frame_path_parts{ii}, '+nosnoc')
            id_complete = true;
        end
        error_id = [normalize_part_for_id(last_frame_path_parts{ii}) ':' error_id];
        ii = ii-1;
    end

    error_struct.identifier = error_id;
    error_struct.stack = stack;
    
    error(error_struct);
end


function [part] = normalize_part_for_id(part)
    part = strip(part, 'left', '+'); % remove leading + if it exists.
    part = char(part); % return char for concatenation.
end
