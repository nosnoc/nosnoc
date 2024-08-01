function model = reformulation(cls_model, opts)
% The main time freezing time freezing reformulation method for Complementarity Lagrangian Systems.
% 
% Args:
%     cls_model (nosnoc.model.Cls): The complementarity Lagrange system to be transformed.
%     opts (nosnoc.Options): Options object.
%
% Returns:
%     nosnoc.model.Heaviside | nosnoc.model.Pss: Returns a model .
    if any(cls_model.e ~= 0)
        model = nosnoc.time_freezing.cls_elastic(cls_model, opts);
    else
        if opts.tf_multicontact
            model = nosnoc.time_freezing.cls_inelastic_multicontact(cls_model, opts);
        else
            model = nosnoc.time_freezing.cls_inelastic(cls_model, opts);
        end
    end
end

