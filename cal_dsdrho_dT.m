function [dsdrho_T,dsdT_rho] = cal_dsdrho_dT(rho,T,sub)
Cv = cal_capacity(rho,T,sub);
dpdT = EOS_PR_dT(rho,T,sub);
dsdrho_T = -dpdT./rho.^2;
dsdT_rho = Cv./T;
end