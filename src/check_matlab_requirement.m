function check_matlab_requirement()
    % 9.11 is R2021b
    % https://www.mathworks.com/support/requirements/previous-releases.html
    v = version();
    c = strsplit(v, '.');
    major_version = str2double(c{1});
    if major_version < 9
        error('Matlab version is incompatible with nosnoc, requires R2021b or later.');
    else
        if major_version == 9
            minor_version = str2double(c{2});
            if minor_version < 11
                error('Matlab version is incompatible with nosnoc, requires R2021b or later.');
            end
        end
    end
end
