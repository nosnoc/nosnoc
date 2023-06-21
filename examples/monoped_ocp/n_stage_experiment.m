import casadi.*

stages = [30,35,40,45,50,55,60,65,70,75,80,85,90];
mode = [0,1];
[X,Y] = meshgrid(stages,mode);
experiments = num2cell([X(:) Y(:)],2)


parfor ii=1:length(experiments)
    experiments(ii)
end
