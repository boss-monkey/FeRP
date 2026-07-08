% Calculate energy
function e = cal_energy(rho,T,sub)

    [a,b,M,dadT] = cal_PR(T,sub);
    [Tc,~,~,~,~,cv_c,n,e_c] = cal_critical_parameter(sub);

    e_ideal = cv_c./((n+1).*Tc.^n).*(T.^(n+1)-Tc.^(n+1));
    
    V = M./rho;
    K0 = 1/sqrt(8)./b.*log((V+(1-sqrt(2)).*b)./(V+(1+sqrt(2)).*b));
    e = e_c+e_ideal + (a - T.*dadT).*K0./M;


    end