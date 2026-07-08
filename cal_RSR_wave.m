function [State_a, State_b] = cal_RSR_wave(W, Pa_guess, Pb_guess,sub)

    tol = 1e-0;       
    maxStep = 500;   
    relax_factor = 1; 
    
    X_k = [Pa_guess; Pb_guess]; %
    
    % --- L1 ---
    for iter = 1:maxStep
        
        pa = X_k(1);
        pb = X_k(2);
        
        % ---  (C1 ) ---
        rhoa = cal_rho_isentropic(W, pa,sub); 

        [ca, ~, ~] = cal_eq_parameter(pa, rhoa,sub); 
        va = 1 / rhoa;
        rhoc_a_sq = (rhoa * ca)^2;
        
        % --- (C2) ---
        Wa = [rhoa;0; pa];
        
        rhob = cal_rho_hugoniot(Wa, pb,sub);
        [cb, ~, ~] = cal_eq_parameter(pb, rhob,sub); 
        vb = 1 / rhob;
        rhoc_b_sq = (rhob * cb)^2;

        j_sq = (pb - pa) / (va - vb);
        
        f1 = j_sq - rhoc_a_sq; % C3 
        f2 = j_sq - rhoc_b_sq; % C4
        F = [f1; f2];
        

        [dcdrho_p_a, dcdp_rho_a, ~, ~, dsdrho_p_a, dsdp_rho_a] = cal_eq_derivative(pa, rhoa,sub);
        [dcdrho_p_b, dcdp_rho_b, dedrho_p_b, dedp_rho_b, ~, ~] = cal_eq_derivative(pb, rhob,sub);
        
      
        Da = -dsdp_rho_a / dsdrho_p_a;
        dFrh_dPb   = -dedp_rho_b + (rhob - rhoa) / (2 * rhoa * rhob);
        dFrh_drhob = -dedrho_p_b + (pa + pb) / (2 * rhob^2);
        
        if abs(dFrh_drhob) < 1e-14
            fprintf('Error: d(Frh)/d(rhob) approaches zero. Implicit derivative cannot be computed. \n');
            J = nan(2,2);
        else
            Db_Pb = -dFrh_dPb / dFrh_drhob;

            % Db_Pa = d(rhob)/d(Pa) | H (C2 对 Pa 的导数)
            % (d(Frh)/d(Pa))_ext + d(Frh)/d(rhob) * Db_Pa = 0
%             de_a_dPa_s = (pa / rhoa^2) * Da; % (来自 d(e)/d(rho)|s * d(rho)/d(P)|s)

%             dFrh_dPa_ext = (1/rhoa - (rhoa+rhob)/(2*rhoa*rhob)) + ... % d(Frh)/d(Pa)
%                            (-pa/rhoa^2 - (pb-pa)/(2*rhoa^2)) * Da + ... % d(Frh)/d(rhoa) * d(rhoa)/d(Pa)
%                            (1) * de_a_dPa_s;                           % d(Frh)/d(ea) * d(ea)/d(Pa)
                       
            dFrh_dPa_ext = (rhob - rhoa)/(2*rhoa*rhob)- (pb-pa)/(2*rhoa^2)* Da;
                          

            Db_Pa = -dFrh_dPa_ext / dFrh_drhob;

            dj_sq_dPa = (-1 / (va - vb)) + ...                         % d(j^2)/d(Pa)
                        (-j_sq / (va - vb)) * (-1/rhoa^2 * Da) + ...  % d(j^2)/d(Va) * d(Va)/d(Pa)
                        (-j_sq / (va - vb)) * (-1/rhob^2) * Db_Pa;    % d(j^2)/d(Vb) * d(Vb)/d(rho) * d(rho)/d(Pa)

            dj_sq_dPb = (1 / (va - vb)) + ...                         % d(j^2)/d(Pb)
                        (-j_sq / (va - vb)) * (-1/rhob^2) * Db_Pb;    % d(j^2)/d(Vb) * d(Vb)/d(rho) * d(rho)/d(Pb)

            d_rhoc_a_sq_dPa = (2 * rhoa^2 * ca * dcdp_rho_a) + ...      % d(rhoc^2)/d(P) | rho
                              (2 * rhoa * ca * (ca + rhoa * dcdrho_p_a)) * Da; % d(rhoc^2)/d(rho) | P * d(rho)/d(P)|s

            d_rhoc_b_sq_dPa = (2 * rhob * cb * (cb + rhob * dcdrho_p_b)) * Db_Pa; % d(rhoc^2)/d(rho) | P * d(rho)/d(Pa)

            d_rhoc_b_sq_dPb = (2 * rhob^2 * cb * dcdp_rho_b) + ...      % d(rhoc^2)/d(P) | rho
                              (2 * rhob * cb * (cb + rhob * dcdrho_p_b)) * Db_Pb; % d(rhoc^2)/d(rho) | P * d(rho)/d(P)|H

            J11 = dj_sq_dPa - d_rhoc_a_sq_dPa;
            J12 = dj_sq_dPb;
            J21 = dj_sq_dPa - d_rhoc_b_sq_dPa;
            J22 = dj_sq_dPb - d_rhoc_b_sq_dPb;
            
            J = [J11, J12; J21, J22];
        end
        
        % 检查 J 是否奇异
        if any(isnan(J), 'all') || abs(det(J)) < 1e-12
            fprintf('Error: The Jacobian matrix is singular or invalid. (det=%.2e) \n', det(J));
            State_a = [pa, rhoa, ca];
            State_b = [pb, rhob, cb];
            return;
        end
        
        delta = J \ F; 
        if X_k(1)*0.5 < delta(1) || X_k(2)*0.5 < delta(2)
             delta1 = 0.5*X_k.*abs(delta)./delta;
        else
            delta1 = delta;
        end        
        
        
        X_k_new = X_k - relax_factor * delta1;

        err = norm(F);
        step_size = norm(X_k_new - X_k) / (norm(X_k) + 1e-9);
        if mod(iter, 20) == 0
        fprintf('iter %d: |F| = %.3e, Step = %.3e, Pa = %.4e, Pb = %.4e\n', ...
                 iter, err, step_size, X_k_new(1), X_k_new(2));
        end
        
        X_k = X_k_new;

        if err < tol && step_size < tol
            fprintf('Dual-sonic shock wave calculation successfully converged');
            rhoa = cal_rho_isentropic(W, X_k(1),sub);
            [ca, ~, ~] = cal_eq_parameter(X_k(1), rhoa,sub);
            Wa = [rhoa; W(2); X_k(1)];
            rhob = cal_rho_hugoniot(Wa, X_k(2),sub);
            [cb, ~, ~] = cal_eq_parameter(X_k(2), rhob,sub);
            
            State_a = [X_k(1), rhoa, ca];
            State_b = [X_k(2), rhob, cb];
            return;
        end
    end
    
    fprintf('Error: R-S wave failed to converge after %d iterations.\n', maxStep);
    State_a = [pa, rhoa, ca];
    State_b = [pb, rhob, cb];
end
