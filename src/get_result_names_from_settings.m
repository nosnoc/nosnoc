function names = get_result_names_from_settings(settings)

names = {"x", "v", "z"};
switch settings.dcs_mode
  case "Stewart"
    names = [names, "theta", "lam", "mu"];
  case "Step"
    names = [names, "alpha", "lambda_n", "lambda_p"];
  case "CLS"
    names = [names, "lambda_normal", "lambda_tangent", "y_gap", "gamma", "beta_conic", "gamma_d", "beta_d", "delta_d", "p_vt", "n_vt", "alpha_vt", "x_left_bp",...
             "Y_gap", "Lambda_normal", "Lambda_tangent", "Gamma", "Gamma_d", "Beta_conic", "Beta_d", "Delta_d", "L_vn", "N_vt", "Alpha_vt"];
end