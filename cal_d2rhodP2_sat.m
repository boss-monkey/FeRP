function d2rhodP2 = cal_d2rhodP2_sat(rho, Tsat, dTdP, d2TdP2,sub)
dpdT = EOS_PR_dT(rho, Tsat,sub);
dpdrho = EOS_PR_drho(rho, Tsat,sub);
drhodP = (1-dpdT*dTdP)/dpdrho;

 [d2pdrho2,d2pdTdrho,d2pdT2] = EOS_PR_d2p(rho,Tsat,sub);


U = (1-dpdT*dTdP);
V = dpdrho;

dUdP = -( (d2pdTdrho.*drhodP + d2pdT2.*dTdP).*dTdP + dpdT.* d2TdP2);
dVdP = d2pdrho2.* drhodP +d2pdTdrho.*dTdP;

d2rhodP2 = (dUdP.*V -U.*dVdP)./V.^2;

end