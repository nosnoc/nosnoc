classdef NosnocDimensions < handle
    properties
        N_stages
        N_finite_elements
        n_x
        n_f
        n_u
        n_q
        n_z
        n_s
        n_sys
        n_c_sys
        n_f_sys
        n_p_global
        n_p_time_var
        n_contacts
        n_tangents
        n_t

        n_theta
        n_lambda
        
        n_alpha
        n_lambda_n
        n_lambda_p

        n_quad
        n_aux
        n_dim_contact

        n_theta_step
        n_beta
    end
end
