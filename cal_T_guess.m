function T0 = cal_T_guess(P,sub)
[Tc,Pc,~,omega] = cal_critical_parameter(sub);

    Pr = P/ Pc;
    k = 5.3727 * (1 + omega);
    discriminant = k^2 - 4 * k * log(Pr);
    Tr_guess = 2 * k ./ (k + discriminant.^0.5);

    T0 = Tr_guess * Tc;
end