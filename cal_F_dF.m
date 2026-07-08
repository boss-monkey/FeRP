function [F,dFdp] = cal_F_dF(W_1, p_star,sub)
%Compute compression branch and its derivatives.
if nargin < 5
    sub = 0;
end

rho = W_1(1,:);
p   = W_1(3,:);

rho_star = cal_rho_hugoniot(W_1, p_star,sub);

j = ((p_star-p).*(rho.*rho_star)./(rho_star-rho+(1e-10))).^0.5;
F =  j.*(1./rho_star-1./rho);

[~,~,dedrho_p,dedp_rho]= cal_eq_derivative(p_star,rho_star,sub);
drho =  -dedp_rho + (rho_star-rho)./(2.*rho.*rho_star);
dp =  -dedrho_p+(p+p_star)./(2.*rho_star.^2);
drhodp = -drho./dp;

dFdp  = -1./2./j-j/2./rho_star.^2*drhodp;


end