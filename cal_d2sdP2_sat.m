function d2sdP2 = cal_d2sdP2_sat(rho,Tsat,dTdP,d2TdP2,sub)
[d2sdrho2,d2sdTdrho,d2sdT2] = cal_d2s(rho,Tsat,sub);
[dsdrho_T,dsdT_rho] = cal_dsdrho_dT(rho,Tsat,sub);
d2rhodP2 = cal_d2rhodP2_sat(rho, Tsat, dTdP, d2TdP2,sub);
drhodP = cal_drhodP_sat(rho, Tsat, dTdP,sub);

d2sdP2 = d2sdrho2.*drhodP.^2 +2*d2sdTdrho.*dTdP.*drhodP +d2sdT2.*dTdP.^2+dsdrho_T.*d2rhodP2+dsdT_rho.*d2TdP2;

end