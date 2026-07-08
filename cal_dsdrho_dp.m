function [dsdrho_p,dsdp_rho] = cal_dsdrho_dp(rho,T,sub)

dpdT = EOS_PR_dT(rho,T,sub);
dpdrho = EOS_PR_drho(rho,T,sub);
Cv = cal_capacity(rho,T,sub);

dsdrho_p = -Cv./T.*dpdrho./dpdT-1./rho.^2.*dpdT;
dsdp_rho = Cv./dpdT./T;
end