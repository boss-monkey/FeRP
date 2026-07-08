function [J,dJdp] = cal_J_dJ_nonclassical(W_1,type,States,p_star,sub)
%Calculate the expansion branch and its derivatives.

    if nargin < 5
        sub = 0;
    end


state_a = States.a;
state_b = States.b;
state_c = States.c; 
state_d = States.d;

pa = state_a(1);
rhoa = state_a(2);
ca = state_a(3);
pb = state_b(1);
rhob = state_b(2);

pc = state_c(1);
pd = state_d(1);

if type == 1
    J = cal_J(W_1,p_star,sub);
    dJdp =  cal_dJdp(W_1, p_star,sub);
else
    if p_star > pd
        J = cal_J(W_1,p_star, sub);
        dJdp =  cal_dJdp(W_1, p_star, sub); %R
    elseif p_star <= pd &&  p_star > pb
        [pa_sol, ~, j_sol] = cal_RS_wave(W_1, p_star, pc, sub);
        J1 = cal_J(W_1,pa_sol, sub);
        J2 = (pa_sol - p_star)./j_sol;
        J = J1+J2;
        dJdp =  cal_dJ_RS(W_1, p_star, pa_sol, sub);  %R-S
    else
        J1 =  cal_J(W_1,pa,sub);
        J2 = (pa-pb)./(rhoa*ca);
        J3 = cal_J([rhob;0;pb],p_star,sub);
        J = J1+J2+J3;
        dJdp = cal_dJdp([rhob;0;pb],p_star,sub);  %R-S-R
    end
end

 
end
