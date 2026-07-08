function rho = cal_rho_isentropic(W_1,p_star,sub,maxStep,tol)

    if nargin < 5
        tol = 1e-3; 
    end
    if nargin < 4
        maxStep = 100;
    end

rho0 = W_1(1,:);
p   = W_1(3,:);

[~,~,s0,~,~,~]= cal_eq_parameter(p,rho0,sub);

rho_current = rho0;


for i = 1:maxStep
    [~,~,s,~,~,~] = cal_eq_parameter(p_star, rho_current, sub);
    [~,~,~,~,dsdrho_p] = cal_eq_derivative(p_star, rho_current, sub);

    F = s - s0;

    if max(abs(F)) < tol
        break;
    end

    min_deriv = 1e-12;
    dsdrho_p(abs(dsdrho_p) < min_deriv) = sign(dsdrho_p(abs(dsdrho_p) < min_deriv)) * min_deriv;

    raw_step = F ./ dsdrho_p;

    max_change = 0.4 * rho_current; 

    safe_step = sign(raw_step) .* min(abs(raw_step), max_change);

    rho_next = rho_current - safe_step;

    neg_mask = rho_next <= 1e-6;
    if any(neg_mask)

        rho_next(neg_mask) = rho_current(neg_mask) * 0.5;
    end

    rho_current = rho_next;
end

rho = rho_current;

if i == maxStep
    warning('[cal_rho_isentropic] Maximum iterations reached. Residual: %e', max(abs(F)));
end   



end