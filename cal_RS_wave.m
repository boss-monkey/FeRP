function [Pa_sol, rhoa_sol, j_sol] = cal_RS_wave(W, Pb, Pa_guess,sub)

    tol = 1;
    maxStep = 100;
    relax_factor = 0.9;
    
    Pa = Pa_guess;
    
     disp(['Iteration into R-S wave, post-wave pressure = ' num2str(Pb) 'pa']); 

    for iter = 1:maxStep
        
        % ---  (C1) ---
        rhoa = cal_rho_isentropic(W, Pa,sub); 
        [ca, ~, ~] = cal_eq_parameter(Pa, rhoa,sub); 
        va = 1 / rhoa;
        ja2 = (rhoa * ca)^2; % C3 

        % ---(C2) ---
        Wa = [rhoa; W(2); Pa];
        rhob = cal_rho_hugoniot(Wa, Pb,sub); 

        vb = 1 / rhob;
        
        if abs(va - vb) < 1e-14
            fprintf('Error: Va and Vb too close. \n');
            Pa_sol = NaN; rhoa_sol = NaN; j_sol = NaN; return;
        end
        
        j2 = (Pb - Pa) / (va - vb); 
        
        F = j2 - ja2; 
        

        %  f = dF/dPa (J_11)
        
        [dcdrho_p_a, dcdp_rho_a, ~, ~, dsdrho_p_a, dsdp_rho_a] = cal_eq_derivative(Pa, rhoa,sub);
        [~, ~, dedrho_p_b, ~, ~, ~] = cal_eq_derivative(Pb, rhob,sub);
        
        % Da = d(rho_a)/d(P_a) | S_L
        Da = -dsdp_rho_a / dsdrho_p_a;
        
        % Db_Pa = d(rhob)/d(Pa) | H
        % d(Frh)/d(Pa)_ext + d(Frh)/d(rhob) * Db_Pa = 0
        dFrh_drhob = -dedrho_p_b + (Pa + Pb) / (2 * rhob^2); 
        
        if abs(dFrh_drhob) < 1e-12
             fprintf('Error: d(Frh)/d(rhob) approaches zero. Implicit derivative cannot be computed.\n');
             Pa_sol = NaN; rhoa_sol = NaN; j_sol = NaN; return;
        end
        
        dedpa_s = (Pa / rhoa^2) * Da; % d(ea)/d(Pa)|s
        
        dFrh_dPa_ext = (1/rhoa - (rhoa+rhob)/(2*rhoa*rhob)) + ... % d(Frh)/d(Pa)
                       (-Pa/rhoa^2 - (Pb-Pa)/(2*rhoa^2)) * Da + ... % d(Frh)/d(rhoa) * d(rhoa)/d(Pa)
                       (1) * dedpa_s;                           % d(Frh)/d(ea) * d(ea)/d(Pa)
        
        Db_Pa = -dFrh_dPa_ext / dFrh_drhob;

        % ---  J_11 = dF/dPa ---
        
        % d(j^2)/d(Pa) 
        dj2_dPa_ext = (-1 / (va - vb)) + ...                         % d(j^2)/d(Pa)
                      (-j2 / (va - vb)) * (-1/rhoa^2 * Da);     % d(j^2)/d(Va) * d(Va)/d(Pa)
        dj2_dPa_imp = (-j2 / (va - vb)) * (-1/rhob^2) * Db_Pa;  % d(j^2)/d(Vb) * d(Vb)/d(rho) * d(rho)/d(Pa)
        dj2_dPa = dj2_dPa_ext + dj2_dPa_imp;

        % d(rho_a*c_a)^2 / d(P_a) 
        dja2_dPa = (2 * rhoa^2 * ca * dcdp_rho_a) + ...      % d(rhoc^2)/d(P) | rho
                          (2 * rhoa * ca * (ca + rhoa * dcdrho_p_a)) * Da; % d(rhoc^2)/d(rho) | P * d(rho)/d(P)|s

        f = dj2_dPa - dja2_dPa;

        if abs(f) < 1e-14
            fprintf('Error: Derivative f approaches zero (Pa=%.3e). \n', Pa);
            Pa_sol = NaN; rhoa_sol = NaN; j_sol = NaN; return;
        end
        

        delta = F / f;
        if abs(delta) > Pa*0.5
             delta1 = Pa*0.5* abs(delta) ./delta;
        else
            delta1 =delta;
        end       
        
        Pa_new = Pa - relax_factor * delta1;
        

        
        err = abs(F);
        step_size = abs(Pa_new - Pa) / (abs(Pa));        
        Pa = Pa_new;

        if err < tol && step_size < tol
            disp(['R-S wave calculation successful, front-wave pressure = ' num2str(Pa) 'pa']); 
            Pa_sol = Pa;
            rhoa_sol = rhoa; 
            j_sol = j2^0.5;
            return;
        end
    end
    
    fprintf('Error: R-S wave failed to converge after %d iterations.\n', maxStep);
    Pa_sol = NaN; rhoa_sol = NaN; j_sol = NaN;
end