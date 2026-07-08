function [dCpdrho, dCpdT] = cal_dCpdrho_dT(rho, T,sub)

[d2pdrho2, d2pdTdrho, d2pdT2] = EOS_PR_d2p(rho, T,sub);
dpdrho = EOS_PR_drho(rho, T,sub); 
dpdT   = EOS_PR_dT(rho, T,sub);   

[dCvdrho,dCvdT] = cal_dCvdrho_dT(rho, T,sub);

rho2 = rho.^2;
dpdT2 = dpdT.^2;
dpdrho_inv = 1./dpdrho; 

T1_T = dCvdT;
T2_T = dpdT2 .* dpdrho_inv ./ rho2;
T3_T = 2.*T .* dpdT .* d2pdT2 .* dpdrho_inv ./ rho2;
T4_T = -T .* dpdT2 .* d2pdTdrho .* (dpdrho_inv.^2) ./ rho2;

dCpdT = T1_T + T2_T + T3_T + T4_T;


rho3 = rho.^3;
T1_rho = dCvdrho;
T2_rho = -2.*T .* dpdT2 .* dpdrho_inv ./ rho3;
T3_rho = 2.*T .* dpdT .* d2pdTdrho .* dpdrho_inv ./ rho2;
T4_rho = -T .* dpdT2 .* d2pdrho2 .* (dpdrho_inv.^2) ./ rho2;

dCpdrho = T1_rho + T2_rho + T3_rho + T4_rho;

end