function [settings] = fill_in_missing_settings(varargin)
if nargin == 1
    settings = varargin{1};
    model = struct([]);
else
    settings = varargin{1};
    model = varargin{2};
end
    


[default_settings] = default_settings_fesd();
fields_settings = fieldnames(settings);
fields_defualt_settings = fieldnames(default_settings);
fields_model = fieldnames(model);

jj = 0;
for ii = 1:length(fields_defualt_settings)
    if ~ismember(fields_defualt_settings{ii},fields_settings)
        eval(['settings.' fields_defualt_settings{ii} '= default_settings.' fields_defualt_settings{ii} ';' ])
        jj = jj+1;
    end
end

fprintf('Added %d missing fields to the user provided settings.\n',jj);


%% The updated settings field.
fields_settings = fieldnames(settings);

jj = 0;
for ii = 1:length(fields_model)
    if ismember(fields_model{ii},fields_settings)
        eval(['settings.' fields_model{ii} '= model.' fields_model{ii} ';' ])
        jj = jj+1;
    end
end
fprintf('Total of %d default settings fields overwritten by model data.\n',jj);



%% Correct contradicting and missing settings
settings = refine_settings(settings);



end

