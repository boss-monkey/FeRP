function [dCvdrho,dCvdT] = cal_dCvdrho_dT(rho, T,sub)

    [~, b, M, ~, d2adT2] = cal_PR(T,sub);
    [Tc,~,~,~,~,cv_c,n] = cal_critical_parameter(sub);
    
    cv_ideal = cv_c.*(T./Tc).^n;
    dCv0dT = n.*cv_ideal./T;
    V = M ./ rho; 
    K0 = 1./(sqrt(8)*b).*log((V+(1-sqrt(2)).*b)./(V+(1+sqrt(2)).*b));
    d3adT3 = d2adT2 .* (-1.5 ./ T);
    d_residual_dT = (K0 ./ M) .* (d2adT2 + T .* d3adT3);
   
    dCvdT = dCv0dT - d_residual_dT;
    
    [~,~,d2pdT2] = EOS_PR_d2p(rho,T,sub);
    dCvdrho = -1./rho.^2.*T.*d2pdT2;

end