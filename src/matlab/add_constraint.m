function [problem] = add_constraint(problem, g, lbg, ubg, idx)
    problem.g = {problem.g{:}, g};
    problem.lbg = [problem.lbg; lbg];
    problem.ubg = [problem.lbg; ubg];
    % This is a bit of a hack, to avoid a bunch of if statements (maybe, currently would only be one)
    if exist('idx','var')
        ind_set = problem.(strcat('ind_', idx));
        problem.(strcat('ind_', idx)) = [ind_set, length(problem.lbg)];
    end
end
