function [problem] = add_variable(problem, w, w0, lbw, ubw, idx)
    %ADD_VARIABLE Adds a primal variable with upper and lower bounds and initialization to the problem.
    dims = [length(w), length(w0), length(lbw), length(ubw)];
    if(~all(dims == dims(1)))
        error("dimension mismatch, with dims: w: %d, w0: %d, lbw: %d, ubw: %d", dims(1), dims(2), dims(3), dims(4));
    end
    problem.w = {problem.w{:}, w};
    % This is a bit of a hack, to avoid a bunch of if statements.
    if exist('idx','var')
        ind_set = problem.(strcat('ind_', idx));
        problem.(strcat('ind_', idx)) = [ind_set, length(problem.w0)+1:length(problem.w0)+length(w)];
    end
    problem.w0 = [problem.w0; w0];
    problem.lbw = [problem.lbw; lbw];
    problem.ubw = [problem.ubw; ubw];
end

