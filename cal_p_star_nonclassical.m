function star = cal_p_star_nonclassical(W_R,W_L,type_L,States_L,maxStep,tol)
%Solve intermediate pressure in Riemann problem via Newton's method
    if nargin < 6
        tol = 1; 
    end
    if nargin < 5
        maxStep = 100;
    end

    disp(' '); 
    disp('————Start calculating the Riemann problem ————')

    p_L = W_L(3,:);
    p_R = W_R(3,:);
    u_R = W_R(2,:);
    u_L = W_L(2,:);

    sub0 = 0;
    sub1 = 1;

    p_star = (p_L .* p_R).^0.5; 
    

    N = length(p_star);
    lambda = 0.5 * ones(1, N); 
    G_old = zeros(1, N);    
    

    converged_mask = false(1, N);

    for i = 1:maxStep
        [L, dL] = cal_J_dJ_nonclassical(W_L, type_L, States_L, p_star, sub0);
        [R, dR] = cal_F_dF(W_R, p_star, sub1);
        
        % G(p) = u_L - u_R + f_L(p) + f_R(p) = 0
        G = (u_L - u_R) + L + R;
        dGdp = dL + dR;
        

        dGdp(abs(dGdp) < 1e-12) = 1e-12;

        delta = G ./ dGdp;
        
        if i > 1
            oscillation_mask = (G .* G_old < 0);
            
            diverge_mask = (abs(G) > abs(G_old));
            
            accelerate_mask = ~(oscillation_mask | diverge_mask);
            
            lambda(oscillation_mask | diverge_mask) = max(0.01, lambda(oscillation_mask | diverge_mask) * 0.5);
            
            lambda(accelerate_mask) = min(1.0, lambda(accelerate_mask) * 1.2);
        end
        

        delta_lim = delta;
        limit_mask = abs(delta) > 0.5 * p_star;
        delta_lim(limit_mask) = 0.5 * p_star(limit_mask) .* sign(delta(limit_mask));

        p_star_new = p_star - lambda .* delta_lim;
        
        

        current_error = abs(G); 
        

        mask_update = ~converged_mask; 
        p_star(mask_update) = p_star_new(mask_update);
        G_old = G; 

        converged_mask = converged_mask | (current_error < tol);
        

        if all(converged_mask)
            disp(['Converged at step ', num2str(i)]);
            break;
        end
    end
    
    if i == maxStep
        disp('Warning: Maximum iterations reached. Results may be inaccurate.');
    end


    [L, ~] = cal_J_dJ_nonclassical(W_L, type_L, States_L, p_star, sub0);
    [R, ~] = cal_F_dF(W_R, p_star, sub1);
    
    star = real([p_star; L; R]);
    

    disp(' '); 
    disp(['————Intermediate pressure (mid-point) = ' num2str(p_star(ceil(end/2))) ' pa————'])
    disp(' '); 
end