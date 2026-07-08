function [Tsat,rho_v,rho_l,V_v, V_l] = cal_T_saturation(P,sub)
    Ru = 8.31443;
    T = cal_T_guess(P,sub);
    iter = 0;
    maxStep = 100; 
    converged = false;
    
    while ~converged && iter < maxStep
         [a,b,M,dadT] = cal_PR(T,sub);
        
        A = a * P / (Ru^2 * T.^2);
        B = b * P / (Ru * T);
        
        coeffs = [1, -(1 - B), (A - 2*B - 3*B.^2), -(A*B - B.^2 - B.^3)];
        Z_roots = roots(coeffs);
        
        Z_real = real(Z_roots(abs(imag(Z_roots)) < 1e-8));
        Z_real = sort(Z_real);
        
        Z_l = Z_real(1);
        Z_v = Z_real(end);
        
        lnphiv = cal_lnphi(A,B,Z_v);
        lnphil = cal_lnphi(A,B,Z_l);
        
        W = lnphiv - lnphil;
        
        dlnphiv_dT = cal_dlnphidT(A,B,Z_v,a,dadT,T);
        dlnphil_dT = cal_dlnphidT(A,B,Z_l,a,dadT,T);
        dWdT = dlnphiv_dT - dlnphil_dT;
        
        
        if abs(dWdT) < 1e-12
            dWdT = sign(dWdT) * 1e-12; 
            if dWdT == 0, dWdT = 1e-12; end
        end

        delta = -W / dWdT;
        max_change = 0.2 * T;
        delta_safe = sign(delta) * min(abs(delta), max_change);
        
        T_new = T + delta_safe;

        if T_new <= 1e-2
            T_new = T * 0.8;
        end
        
        if abs(delta_safe) < 1e-6 * T && abs(W) < 1e-8
            converged = true;
        end
        
        T = T_new;
        iter = iter + 1;
    end
    
    V_l = Z_l * Ru * T ./ P;
    V_v = Z_v * Ru * T ./ P;
    rho_l = M./V_l;
    rho_v = M./V_v;
    Tsat = T;
    
    if iter == maxStep
        warning('Saturation T calculation did not converge.');
    end
end