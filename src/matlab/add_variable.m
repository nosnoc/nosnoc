function [problem] = add_variable(problem, w, w0, lbw, ubw, idx)

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

