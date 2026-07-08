function [dcdrho_T,dcdT_rho] = cal_dcdrho_dT(rho,T,sub)


    Cv =  cal_capacity(rho,T,sub);
    dpdrho = EOS_PR_drho(rho,T,sub);
    dpdT = EOS_PR_dT(rho,T,sub);

    c = (dpdrho + T./Cv./rho.^2.*dpdT.^2).^0.5;
    [d2pdrho2,d2pdTdrho,d2pdT2] = EOS_PR_d2p(rho,T,sub);
    
    [dCvdrho,dCvdT] = cal_dCvdrho_dT(rho, T,sub);
    dcdrho_T = 1/2./c.*(d2pdrho2 - 2*T./Cv./rho.^3 .*dpdT.^2 + 2*T./Cv./rho.^2.* dpdT.*d2pdTdrho - T./Cv.^2./rho.^2 .*dCvdrho.*dpdT.^2 );
    dcdT_rho = 1/2./c.*(d2pdTdrho + 1./Cv./rho.^2.*dpdT.^2 + 2*T./Cv./rho.^2.* dpdT.*d2pdT2- T./Cv.^2./rho.^2 .*dCvdT.*dpdT.^2 );
    

end