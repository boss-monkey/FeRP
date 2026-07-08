function dJ = cal_dJ_RS(W, Pb, Pa, sub)

        
        % --- (C1 ) ---
        rhoa = cal_rho_isentropic(W, Pa, sub); 
        [ca, ~, ~] = cal_eq_parameter(Pa, rhoa, sub); 
        va = 1 / rhoa;
        ja = rhoa * ca;

        % ---  (C2 ) ---
        Wa = [rhoa; 0; Pa];
        rhob = cal_rho_hugoniot(Wa, Pb, sub); 

        vb = 1 / rhob;
        j2 = (Pb - Pa) / (va - vb); 

        [dcdrho_p_a, dcdp_rho_a, ~, ~, dsdrho_p_a, dsdp_rho_a] = cal_eq_derivative(Pa, rhoa, sub);
        [~, ~, dedrho_p_b, dedp_rho_b, ~, ~] = cal_eq_derivative(Pb, rhob, sub);

        % Da = d(rho_a)/d(P_a) | S_L
        Da = -dsdp_rho_a / dsdrho_p_a;
        
        % Db_Pa = d(rhob)/d(Pa) | H
        % d(Frh)/d(Pa)_ext + d(Frh)/d(rhob) * Db_Pa = 0
        dFrh_dp_b = -dedp_rho_b + (rhob - rhoa) / (2 * rhob*rhoa);         
        dFrh_drhob = -dedrho_p_b + (Pa + Pb) / (2 * rhob^2); 
        drhodp_b = -dFrh_dp_b/dFrh_drhob;

        
        dedpa_s = (Pa / rhoa^2) * Da; % d(ea)/d(Pa)|s
        
        dFrh_dPa_ext = (1/rhoa - (rhoa+rhob)/(2*rhoa*rhob)) + ... % d(Frh)/d(Pa)
                       (-Pa/rhoa^2 - (Pb-Pa)/(2*rhoa^2)) * Da + ... % d(Frh)/d(rhoa) * d(rhoa)/d(Pa)
                       (1) * dedpa_s;                           % d(Frh)/d(ea) * d(ea)/d(Pa)
        
        Db_Pa = -dFrh_dPa_ext / dFrh_drhob;

        % --- 3. 计算 J_11 = dF/dPa ---
        
        % d(j^2)/d(Pa) (全导数)
        dj2dPa_ext = (-1 / (va - vb)) + ...                         % d(j^2)/d(Pa)
                      (-j2 / (va - vb)) * (-1/rhoa^2 * Da);     % d(j^2)/d(Va) * d(Va)/d(Pa)
        dj2dPa_imp = (-j2 / (va - vb)) * (-1/rhob^2) * Db_Pa;  % d(j^2)/d(Vb) * d(Vb)/d(rho) * d(rho)/d(Pa)
        dj2dPa = dj2dPa_ext + dj2dPa_imp;

        % d(rho_a*c_a)^2 / d(P_a)
        dja2dPa = (2 * rhoa^2 * ca * dcdp_rho_a) + ...      % d(rhoc^2)/d(P) | rho
                          (2 * rhoa * ca * (ca + rhoa * dcdrho_p_a)) * Da; % d(rhoc^2)/d(rho) | P * d(rho)/d(P)|s
                          
        f1 = dj2dPa - dja2dPa;
        f2 = 1 / (va - vb) *(1- j2/rhob^2* drhodp_b);
        dPadPb = f2./f1 ;
        djdP_a = dja2dPa./2/ja;
        dJ = -1./ ja *(1+(Pa-Pb)/ja*djdP_a*dPadPb);

end