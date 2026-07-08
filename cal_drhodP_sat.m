function drhodP_sat = cal_drhodP_sat(rho, Tsat, dTdP,sub)
dpdT = EOS_PR_dT(rho, Tsat,sub);
dpdrho = EOS_PR_drho(rho, Tsat,sub);
drhodP_sat = (1-dpdT*dTdP)/dpdrho;
end