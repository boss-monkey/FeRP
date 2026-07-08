function [dedrho_p,dedp_rho] = cal_dedrho_dp(rho,T,sub)

p = EOS_PR(rho,T,sub);
dpdT = EOS_PR_dT(rho,T,sub);
dpdrho = EOS_PR_drho(rho,T,sub);
Cv = cal_capacity(rho,T,sub);

dedrho_p = -Cv.*dpdrho./dpdT-1./rho.^2.*(T.*dpdT-p);
dedp_rho = Cv./dpdT;
end