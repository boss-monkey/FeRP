% Calculate specific heat capacity
function cv = cal_capacity(rho,T,sub)

    [~,b,M,~,d2adT2] = cal_PR(T,sub);
    [Tc,~,~,~,~,cv_c,n] = cal_critical_parameter(sub);
    
    cv_ideal = cv_c.*(T./Tc).^n;

    V = M./rho;
    K0 = 1./(sqrt(8)*b).*log((V+(1-sqrt(2)).*b)./(V+(1+sqrt(2)).*b));
    cv = cv_ideal - K0./M.*T.*d2adT2;
end