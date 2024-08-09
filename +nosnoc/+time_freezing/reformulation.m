function model = reformulation(cls_model, opts)
    if any(cls_model.e ~= 0)
        model = nosnoc.time_freezing.cls_elastic(cls_model, opts);
    else
        if opts.time_freezing_Heaviside_lifting
            if opts.dcs_mode == DcsMode.Stewart
                warning('nosnoc: time_freezing:reformulation:wrong_dcs_mode',...
                    'Dcs mode is set to Stewart but using multicontact timefreezing, setting it to Heaviside.')
                opts.dcs_mode = DcsMode.Heaviside;
            end
            model = nosnoc.time_freezing.cls_inelastic_multicontact(cls_model, opts);
        else
            model = nosnoc.time_freezing.cls_inelastic(cls_model, opts);
        end
    end
end

