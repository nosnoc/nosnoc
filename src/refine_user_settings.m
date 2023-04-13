function [settings] = refine_user_settings(settings)
settings_bkp = settings;
unfold_struct(settings,'caller')
settings = settings_bkp;

end

