%This program is figure 4 in the paper.
function star = cal_p_star_nonclassical(W_R,W_L,type_L,States_L,maxStep,tol)
disp(' '); 
disp('————Start calculating the Riemann problem————')

if nargin < 7
    maxStep = 100;
    tol = 5;
end
p_L = W_L(3,:);
p_R = W_R(3,:);
p_star =(p_L.*p_R).^0.5; %Initial value
   
u_R   = W_R(2,:);
u_L   = W_L(2,:);


for i = 1:maxStep

    [L,dL] = cal_J_dJ_nonclassical(W_L,type_L,States_L, p_star);
    [R,dR] = cal_F_dF(W_R, p_star);

%%    
    G = (u_L-u_R)+L+R;
    dGdp = dL+dR;
    delta = G./(dGdp);
    
    if abs(delta)>0.5*p_star
        delta1 = 0.5*p_star*delta/abs(delta);
    else
        delta1 = delta;
    end
    p_star_new = p_star - delta1*0.4;
    
    p_star  = p_star_new;
   if max(abs(delta)) < tol
        break;
   end

end
star = real([p_star;L;R]);
disp(' '); 
disp(['————Intermediate pressure = ' num2str(p_star) 'pa, calculation successful————'])
disp(' '); 
end
