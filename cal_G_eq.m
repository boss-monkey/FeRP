function G = cal_G_eq(p,rho)

 [dcdrho_p,dcdp_rho] = cal_eq_derivative(p,rho);
 c = cal_eq_parameter(p,rho);
 dcdrho_s = c.^2.*dcdp_rho +dcdrho_p;
 
 G = 1+rho./c.*dcdrho_s;

end