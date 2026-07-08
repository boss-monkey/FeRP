function dedrho_T = cal_dedrho_T(rho,T,sub)
p = EOS_PR(rho,T,sub);
dpdT = EOS_PR_dT(rho,T,sub);
dedrho_T =(p-T.*dpdT)./rho.^2;
end