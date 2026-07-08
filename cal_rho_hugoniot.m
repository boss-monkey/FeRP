function rho = cal_rho_hugoniot(W_1,p_star,sub)

    if nargin < 7
        tol = 1e-4; 
    end
    if nargin < 6
        maxStep = 100;
    end


    rho0 = W_1(1,:); 
    p0   = W_1(3,:); 
     [~, e] = cal_eq_parameter(p0, rho0, sub);

    rho_current = rho0;

    for i = 1:maxStep
        [~, e_star] = cal_eq_parameter(p_star, rho_current, sub);
        [~, ~, dedrho_p] = cal_eq_derivative(p_star, rho_current, sub);
        F = e - e_star + p0./rho0 - p_star./rho_current + ...
            (p_star - p0) .* (rho0 + rho_current) ./ (2 * rho0 .* rho_current);
        f = -dedrho_p + (p0 + p_star) ./ (2 .* rho_current.^2);

        min_deriv = 1e-12;
        f(abs(f) < min_deriv) = sign(f(abs(f) < min_deriv)) * min_deriv;

        raw_delta = F ./ f;

        max_change = 0.4 * rho_current;

        safe_delta = sign(raw_delta) .* min(abs(raw_delta), max_change);

        rho_new = rho_current - safe_delta;

        neg_mask = rho_new <= 1e-6;
        if any(neg_mask)
            rho_new(neg_mask) = rho_current(neg_mask) * 0.5;
            safe_delta(neg_mask) = rho_current(neg_mask) * 0.5; 
        end
        

        if max(abs(safe_delta)) < tol
             rho_current = rho_new; 
             break;
        end

        rho_current = rho_new; 
    end
    
    rho = real(rho_current);
    
    if i == maxStep
        warning('[cal_rho_hugoniot] Max step reached. Max Residual F: %e', max(abs(F)));
    end   
  
end

