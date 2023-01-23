function print_casadi_vector(x, fileID)
    if ~exist('fileID')
        fileID = 1;
    end
    for i = 1:length(x)
        fprintf(fileID, strcat(formattedDisplayText(x(i)), "\n"));
    end
end
