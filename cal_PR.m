%Calculate the coefficients of the PR equation and the temperature-related derivatives

function [a,b,M,dadT,d2adT2,c1,c2,c3] = cal_PR(T,sub)
%%
    Ru = 8.31443;
    
% Dodecane
    
    [Tc,pc,~,omega,M] = cal_critical_parameter(sub);
    
    c =  0.37464 + 1.54226*omega - 0.26992*omega^2;
    b = 0.077796073903888*Ru.*Tc./pc;
    a = 0.457235528921382*(Ru.*Tc).^2/pc.*(1+c*(1-(T/Tc).^0.5)).^2;

%%

    G = c*(T/Tc).^0.5./(1+c*(1-(T/Tc).^0.5));    
    dadT = -1./T.*(a.*G);
    d2adT2 = 0.457235528921382*Ru.^2./(2*T).* ...
             (c*(1+c)*Tc/pc.*(Tc./T).^0.5);    



     c1 =  0.457235528921382*Ru.^2.* ...
            (Tc.^2./pc.*c.^2./Tc);
     c2 =  0.457235528921382*Ru.^2.* ...
            (Tc.^2./pc.*  2.*(c+c.^2)./ ((Tc).^0.5));         
     c3 =  0.457235528921382*Ru.^2.* ...
             (Tc.^2./pc.*  (1+c).^2);

end

