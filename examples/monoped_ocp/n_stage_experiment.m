function n_stage_experiment(name, stages)
    import casadi.*
    
    mode = [0,1];
    [X,Y] = meshgrid(stages,mode);
    experiments = num2cell([X(:) Y(:)],2)
    nexp = length(experiments);
    results = cell(nexp, 0);
    stats = cell(nexp, 0);

    addAttachedFiles(gcp,["robot_model_files/robot_model_kinematics.m"])
    parfor ii=1:length(experiments)
        experiment = experiments{ii};
        stages = experiment(1);
        mode = experiment(2);
        if mode
            [result,stat] = monoped_model(stages,false,false);
        else
            [result,stat] = monoped_stewart_model(stages,false,false);
        end
        results{ii} = result;
        stats{ii} = stat;
    end

    filename = sprintf([name,'_%s'], datetime('now','Format','yyyyMMddHHmmss'));

    output.experiments = experiments;
    output.results = results;
    output.stats = stats;

    save(filename, "output");
end

