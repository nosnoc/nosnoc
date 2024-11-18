function warning(warn_id, msg)
% nosnoc specific wrapper for warning which builds the full warn_id automagically.
% This may break if for some reason it is called outside of `+nosnoc` or its subdirectories
% or if the user has a +nosnoc directory on the path somehow (Though this is unlikely due
% to how we walk down the path).
%
% TODO:
%       @anton Why does matlab not allow passing a stack to warning. As such I am not sure this is the best way of doing this.
%       We may just want a method that gets the warn_id and use normal warning. This doesn't match the very nice idiomatic use
%       of nosnoc.error however and that is annoying.
    
    stack = dbstack('-completenames', 1);

    [last_frame_path,last_frame_fname,~] = fileparts(stack(1).file);

    last_frame_path_parts = split(last_frame_path, filesep);

    id_complete = false;
    warn_id = [normalize_part_for_id(last_frame_fname) ':' normalize_part_for_id(warn_id)];
    ii = length(last_frame_path_parts);
    while ~id_complete
        if strcmp(last_frame_path_parts{ii}, '+nosnoc')
            id_complete = true;
        end
        warn_id = [normalize_part_for_id(last_frame_path_parts{ii}) ':' warn_id];
        ii = ii-1;
    end
    
    warning(warn_id, ['nosnoc: ' char(msg)]);
end


function [part] = normalize_part_for_id(part)
    part = strip(part, 'left', '+'); % remove leading + if it exists.
    part = char(part); % return char for concatenation.
end
