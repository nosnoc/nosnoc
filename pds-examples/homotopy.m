function [success,stats] = homotopy(prob,sigma_0,comp_tol,slope)
    if ~exist('slope')
        slope = 0.1;
    end
    if ~exist('comp_tol') 
        comp_tol = 1e-9;
    end
    if ~exist('sigma_0') 
        sigma_0 = 1;
    end
    sigma_k = sigma_0;
    all_stats = [];
    while sigma_k >= comp_tol
        prob.p.sigma(1).init = sigma_k;
        stats = prob.solve();
        prob.w.init = prob.w.res;
        all_stats = [all_stats, stats];
        comp_res = full(prob.comp_res_fun(prob.w.res,prob.p.init));
        if stats.return_status == "Search_Direction_Becomes_Too_Small"
            stats.success = 1;
        end
        if ~stats.success
            fprintf(['sigma_k=' num2str(sigma_k) ' comp_res=' num2str(comp_res) '\n']);
        end
        sigma_k = slope*sigma_k;
        if comp_res < comp_tol && stats.success
            break
        end
        if sigma_k == 0
            break
        end
    end

    success = stats.success;
    if ~success
        fprintf([stats.return_status '\n']);
    end
end
