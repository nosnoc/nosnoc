function [problem] = add_variable(problem, w, w0, lbw, ubw, type)

    problem.w = {problem.w{:}, w};
    % This is a bit of a hack, to avoid a bunch of if statements.
    ind_set = problem.(strcat('ind_', type))
    problem.(strcat('ind_', type)) = [ind_set, length(problem.w0)+1:length(problem.w0)+length(w)]
    problem.w0 = [problem.w0; w0];
    problem.lbw = [problem.lbw; w0];
    problem.ubw = [problem.ubw; w0];
end

