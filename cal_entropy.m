% Calculate entropy
function s = cal_entropy(rho,T,sub)

    [~,b,M,dadT] = cal_PR(T,sub);
    [Tc,~,~,~,~,cv_c,n,~,s_c] = cal_critical_parameter(sub);
    Ru = 8.31443;
    R = Ru./M;
    
    s_ideal =  cv_c./(n.*Tc.^n).*(T.^(n)-Tc.^(n));
    V = M./rho;
    K0 = 1/sqrt(8)./b.*log((V+(1-sqrt(2)).*b)./(V+(1+sqrt(2)).*b));
    s = s_c + s_ideal + R.*log((V-b)./M) - K0./M.*dadT;
    
   
end

