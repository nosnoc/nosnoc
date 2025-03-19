function psi_fun = get_psi_fun(relaxation_type, normalize)
    import casadi.*
    import nosnoc.reg_homotopy.*
    a = SX.sym('a',1);
    b = SX.sym('b',1);
    sigma = SX.sym('sigma',1);
    switch relaxation_type
      case {MpccMethod.SCHOLTES_EQ, MpccMethod.SCHOLTES_INEQ}
        psi_mpcc = a.*b-sigma;
        norm = sigma;
      case {MpccMethod.FISCHER_BURMEISTER_EQ,MpccMethod.FISCHER_BURMEISTER_INEQ}
        if normalize
            normalized_sigma = sqrt(2*sigma);
        else
            normalized_sigma = sigma;
        end
        psi_mpcc = a+b-sqrt(a^2+b^2+normalized_sigma^2);
        
      case {MpccMethod.NATURAL_RESIDUAL_EQ,MpccMethod.NATURAL_RESIDUAL_INEQ}
        if normalize
            normalized_sigma = sqrt(4*sigma);
        else
            normalized_sigma = sigma;
        end
        psi_mpcc = 0.5*(a+b-sqrt((a-b)^2+normalized_sigma^2));
      case {MpccMethod.CHEN_CHEN_KANZOW_EQ, MpccMethod.CHEN_CHEN_KANZOW_INEQ}
        alpha = 0.5;
        if normalize
            psi_mpcc = alpha*(a+b-sqrt(a^2+b^2+2*sigma))+(1-alpha)*(a*b-sigma);
        else
            psi_mpcc = alpha*(a+b-sqrt(a^2+b^2+sigma^2))+(1-alpha)*(a*b-sigma);
        end
      case {MpccMethod.STEFFENSEN_ULBRICH_EQ, MpccMethod.STEFFENSEN_ULBRICH_INEQ}
        if normalize
            normalized_sigma = 2/((2/pi)*sin(3*pi/2)+1)*sqrt(sigma);
        else
            normalized_sigma = sigma;
        end
        x = a-b;
        z = x/normalized_sigma;
        y_sin = normalized_sigma*((2/pi)*sin(z*pi/2+3*pi/2)+1);
        psi_mpcc = a+b-if_else(abs(x)>=normalized_sigma,abs(x),y_sin);
      case {MpccMethod.STEFFENSEN_ULBRICH_POLY_EQ, MpccMethod.STEFFENSEN_ULBRICH_POLY_INEQ}
        if normalize
            normalized_sigma = (16/3)*sqrt(sigma);
        else
            normalized_sigma = sigma;
        end
        x = a-b;
        z = x/normalized_sigma;
        y_pol = normalized_sigma*(1/8*(-z^4+6*z^2+3));
        psi_mpcc = a+b- if_else(abs(x)>=normalized_sigma,abs(x),y_pol);
      case {MpccMethod.KANZOW_SCHWARTZ_EQ, MpccMethod.KANZOW_SCHWARTZ_INEQ}
        if normalize
            normalized_sigma = sqrt(sigma);
        else
            normalized_sigma = sigma;
        end 
        a1 = a-normalized_sigma;
        b1 = b-normalized_sigma;
        psi_mpcc = if_else((a1+b1)>=0,a1*b1,-0.5*(a1^2+b1^2));
      case MpccMethod.LIN_FUKUSHIMA
        if normalize
            normalized_sigma = sqrt(sigma);
        else
            normalized_sigma = sigma;
        end
        psi_mpcc1 = a*b-normalized_sigma^2;
        psi_mpcc2 = ((a+normalized_sigma)*(b+normalized_sigma)-normalized_sigma^2);
        psi_mpcc = vertcat(psi_mpcc1, psi_mpcc2);
      case MpccMethod.KADRANI
        if normalize
            normalized_sigma = sqrt(sigma);
        else
            normalized_sigma = sigma;
        end
        psi_mpcc1 = (a-normalized_sigma)*(b-normalized_sigma);
        psi_mpcc2 = -normalized_sigma - a;
        psi_mpcc3 = -normalized_sigma - b;
        psi_mpcc = vertcat(psi_mpcc1, psi_mpcc2, psi_mpcc3);
    end

    psi_fun = Function('psi_fun',{a,b,sigma},{psi_mpcc});
end
