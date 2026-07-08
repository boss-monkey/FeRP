%Caculate the speed of sound
function c = cal_c_rho_T(rho,T,sub)

    Cv = cal_capacity(rho,T,sub);
    dpdrho = EOS_PR_drho(rho,T,sub);
    dpdT = EOS_PR_dT(rho,T,sub);
    c = (dpdrho + T./Cv./rho.^2.*dpdT.^2).^0.5;


end