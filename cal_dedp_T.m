function dedp_T = cal_dedp_T(rho,T)
 p = EOS_PR(rho,T);
dpdT = EOS_PR_dT(rho,T);
dpdrho = EOS_PR_drho(rho,T);

dedp_T = (p-T.*dpdT)./dpdrho./rho.^2;
end