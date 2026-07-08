%This program calculates the shock wave speed.
function U = cal_shock_speed(W_1, m, p_star,sub)

rho = W_1(1,:);
u =  W_1(2,:);
p   = W_1(3,:);
rho_star = cal_rho_hugoniot(W_1, p_star,sub);
j = ((p_star-p).*(rho.*rho_star)./(rho_star-rho-(1e-20))).^0.5;

S = u-m*j./rho;
U = S;
end