function dJdp = cal_dJdp(W_1, p_star,sub)

rho_star = cal_rho_isentropic(W_1, p_star,sub);
c_star= cal_eq_parameter(p_star,rho_star,sub);

dJdp = -1./(rho_star.*c_star);

end