function dFdp = cal_dFdp(W_1, p_star)

rho = W_1(1,:);
p   = W_1(3,:);


rho_star = cal_rho_hugoniot(W_1, p_star);

j = ((p_star-p).*(rho.*rho_star)./(rho_star-rho-(1e-20))).^0.5;
% djdp = j./2./(p_star-p-(1e-20));
% djdrho = -j.*rho./2./rho_star./(rho_star-rho-(1e-20));

[~,~,dedrho_p,dedp_rho]= cal_eq_derivative(p_star,rho_star);
drho =  -dedp_rho + (rho_star-rho)./(2.*rho.*rho_star);
dp =  -dedrho_p+(p+p_star)./(2.*rho_star.^2);
drhodp = -drho./dp;

% dFdp = (1./rho_star-1./rho).*djdp +(- j./(rho_star.^2)+(1./rho_star-1./rho).*djdrho).*drhodp;
dFdp  = -1./2./j-j/2./rho_star.^2*drhodp;


end