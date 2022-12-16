function [problem] = add_constraint(problem, g, lbg, ubg, idx)
    %ADD_CONSTRAINT Adds a general constraint with upper and lower bounds to the problem.

    dims = [length(g), length(lbg), length(ubg)];
    if(~all(dims == dims(1)))
        error("dimension mismatch, with dims: g: %d, lbg: %d, ubg: %d", dims(1), dims(2), dims(3));
    end
    problem.g = {problem.g{:}, g};
    problem.lbg = [problem.lbg; lbg];
    problem.ubg = [problem.ubg; ubg];
    % This is a bit of a hack, to avoid a bunch of if statements (maybe, currently would only be one)
    if exist('idx','var')
        ind_set = problem.(strcat('ind_', idx));
        problem.(strcat('ind_', idx)) = [ind_set, length(problem.lbg)];
    end
end
