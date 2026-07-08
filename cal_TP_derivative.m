function [drhodT_p, drhodP_t, dedT_p] = cal_TP_derivative(p, rho)

    [~, pc, ~, ~, ~] = cal_critical_parameter();
    N = length(p);
    
    drhodT_p = zeros(N, 1);
    drhodP_t = zeros(N, 1);
    dedT_p   = zeros(N, 1);
    

    [~, ~, dedrho_p, dedp_rho, ~, ~] = cal_eq_derivative(p, rho);
    
    T = cal_T_EOS_PR(p, rho);
    
    for i = 1:N
        is_two_phase = false;
        if p(i) < pc * 0.999
            [Tsat, rho_v, rho_l, ~, ~] = cal_T_saturation(p(i));
            if rho(i) > rho_v && rho(i) < rho_l
                is_two_phase = true; 
                drhodT_p(i) = inf;
                drhodP_t(i) = inf;
                dedT_p(i) = inf;
                
            end
        end
        
        if ~is_two_phase

            dpdT_rho = EOS_PR_dT(rho(i), T(i));   
            dpdrho_T = EOS_PR_drho(rho(i), T(i)); 
            drhodT_p(i) = - dpdT_rho / dpdrho_T;
            drhodP_t(i) = 1 / dpdrho_T;
            dedT_p(i) = dedrho_p(i) * drhodT_p(i);
        end
    end
end