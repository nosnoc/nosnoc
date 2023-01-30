function [settings] = fill_in_missing_settings(varargin)
if nargin == 1
    settings = varargin{1};
    model = struct([]);
else
    settings = varargin{1};
    model = varargin{2};
end
    
[default_settings] = default_settings_nosnoc();
fields_settings = fieldnames(settings);
fields_default_settings = fieldnames(default_settings);
fields_model = fieldnames(model);

jj = 0;
for ii = 1:length(fields_default_settings)
    if ~ismember(fields_default_settings{ii},fields_settings)
        eval(['settings.' fields_default_settings{ii} '= default_settings.' fields_default_settings{ii} ';' ])
        jj = jj+1;
    end
end

%% The updated settings field.
fields_settings = fieldnames(settings);

jj = 0;
for ii = 1:length(fields_model)
    if ismember(fields_model{ii},fields_settings)
        eval(['settings.' fields_model{ii} '= model.' fields_model{ii} ';' ])
        jj = jj+1;
    end
end

% %% Correct contradicting and missing settings
settings = refine_settings(settings);

end

