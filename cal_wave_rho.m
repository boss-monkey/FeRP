%This program calculates the flow variables within rarefaction waves.
function w = cal_wave_rho(W_1,m, x_ref, sub, maxStep,tol)
  
    if nargin < 5
        maxStep = 300;
        tol = 1e-3;
    end
    rho0 = W_1(1,:);
    u = W_1(2,:);
    
    N = length(rho0);
    lambda = 0.5 * ones(1, N);
    F_old = inf(1, N);

for j = 1:maxStep
    p0 = cal_p_isentropic(W_1,rho0, sub);
    c0 = cal_eq_parameter(p0,rho0, sub);
    J =  cal_J(W_1,p0, sub);
    [dcdrho,dcdp] = cal_eq_derivative(p0,rho0, sub);
    F = m*u+J-m*x_ref-c0;
    f =  -dcdrho - dcdp.*c0.^2 - c0./rho0;
    delta =  F ./ f;
    
    if j > 1
        mask_osc = (F .* F_old < 0);
        mask_div = (abs(F) >= abs(F_old));
        mask_acc = ~(mask_osc | mask_div);
        
        lambda(mask_osc | mask_div) = max(0.1, lambda(mask_osc | mask_div) * 0.5);
        lambda(mask_acc) = min(1.0, lambda(mask_acc) * 1.2);
    end
    F_old = F;
    
    delta1 = delta;
    mask_limit = abs(delta) > 0.5 * rho0;
    if any(mask_limit)
        delta1(mask_limit) = 0.5 * rho0(mask_limit) .* sign(delta(mask_limit));
    end
    
    rho_new = rho0 - lambda .* delta1;
   
    rho0 = rho_new;

   if max(abs(delta)) < tol
        break;
   end
    
end
   rho1 = rho0;
   p1 =  cal_p_isentropic(W_1,rho1, sub);
   c1 = cal_eq_parameter(p1,rho1, sub);
   u1 = m*c1+x_ref;
   
   w = [rho1;u1;p1];
  disp(['The interface is located within the rarefaction wave, the absolute error is ' num2str(max(abs(delta))) ', max error' num2str(tol)]); 
if j == maxStep
    warning('[cal_wave_rho] Convergence tolerance not met. Results may be inaccurate.');
end   
end
