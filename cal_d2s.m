function [d2sdrho2,d2sdTdrho,d2sdT2] = cal_d2s(rho,T,sub)

cv = cal_capacity(rho,T,sub);
[~,dCvdT] = cal_dCvdrho_dT(rho, T,sub);
[~,d2pdTdrho,d2pdT2] = EOS_PR_d2p(rho,T,sub);
dpdT = EOS_PR_dT(rho,T,sub);

d2sdrho2 = 2./rho.^3 .*dpdT - 1./rho.^2 .*d2pdTdrho;
d2sdTdrho = - 1./rho.^2 .*d2pdT2;
d2sdT2 = dCvdT./T - cv./T.^2;

end