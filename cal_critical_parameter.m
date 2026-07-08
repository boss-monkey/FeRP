function [Tc,pc,rhoc,omega,M,cv_c,n,e_c,s_c] = cal_critical_parameter(sub)
%Coefficients of the Peng-Robinson equation of state
if sub ==0
% % dodecane
     Tc = 658.1;       
    pc = 1.817e6;       
    rhoc = 186;
    omega = 0.574;      
    M =  170.33e-3;

    cv_c = 2.970123153445547e+03;
    n = 0.612914877770408;
    e_c = 6.948178368232952e+05;
    s_c = 1.400759326735013e+03;
else
% % N2
    M =  28.0134e-3;
    Tc = 126.2;
    pc  = 3.40e6;
    rhoc = 313;
    omega  = 0.0372;
    cv_c =7.435175725566218e+02;
    n = 0.086494863530892;
    e_c = 9.308510655790280e+04;
    s_c = 6.236944454704188e+03;
end

  
end