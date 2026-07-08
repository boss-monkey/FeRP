% Calculate specific heat capacity
function cp = cal_capacity_p(rho,T,sub)
    cv = cal_capacity(rho,T,sub);
    dpdrho = EOS_PR_drho(rho,T,sub);
    dpdT = EOS_PR_dT(rho,T,sub);
    cp = cv + T./rho.^2 .*dpdT.^2./dpdrho;
end