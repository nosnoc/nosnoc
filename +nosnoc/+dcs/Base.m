classdef Base < matlab.mixin.Scalar & handle
    properties
        f_x_fun
        f_q_fun
        g_z_fun
        g_alg_fun
        g_path_fun
        g_comp_path_fun
        g_terminal_fun
        f_q_T_fun
        f_lsq_x_fun
        f_lsq_u_fun
        f_lsq_T_fun
    end

    methods(Abstract)
        generate_variables(obj, opts)
        generate_equations(obj, opts)
    end
end
