function p = cal_p_isentropic(W_1,rho1,sub,maxStep,tol)
    
    if nargin < 5
        tol = 1; 
    end
    if nargin < 4
        maxStep = 100;
    end

    rho0 = W_1(1,:);
    p0   = W_1(3,:);
    
    [~,~,s0,~,~,~]= cal_eq_parameter(p0,rho0,sub);

    p_current = p0;
    
    for i = 1:maxStep

        [~,~,s,~,~,~] = cal_eq_parameter(p_current, rho1, sub);
        [~,~,~,~,~,dsdp_rho] = cal_eq_derivative(p_current, rho1, sub);
        F = s0 - s;
        f = -dsdp_rho;

        min_deriv = 1e-12;
        f(abs(f) < min_deriv) = sign(f(abs(f) < min_deriv)) * min_deriv;

        raw_delta = F ./ f;
        

        max_change = 0.5 * p_current;

        safe_delta = sign(raw_delta) .* min(abs(raw_delta), max_change);

        p_new = p_current - safe_delta;

        neg_mask = p_new <= 1e-6;
        if any(neg_mask)
            p_new(neg_mask) = p_current(neg_mask) * 0.5;
            safe_delta(neg_mask) = p_current(neg_mask) * 0.5;
        end

        if max(abs(safe_delta)) < tol
            p_current = p_new;
            break;
        end
        
        p_current = p_new;
    end
    
     p = p_current;
    
    if i == maxStep
        warning('[cal_p_isentropic] Convergence tolerance not met. Max Residual s_err: %e', max(abs(F)));
    end

end