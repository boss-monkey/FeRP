function dJdp = cal_dJdp_nonclassical(W_1,type,States, p_star)


state_a = States.a;
state_b = States.b;
state_c = States.c; 
state_d = States.d;


pb = state_b(1);
rhob = state_b(2);

pc = state_c(1);
pd = state_d(1);

if type == 1
    dJdp =  cal_dJdp(W_1, p_star);
else
    if p_star > pd
        dJdp =  cal_dJdp(W_1, p_star);   %R
    elseif p_star < pd &&  p_star > pb
        [~, ~, j_sol] = cal_RS_wave(W_1, p_star, pc);   %R-S
        dJdp =  -1./j_sol;  
    else
        dJdp = cal_dJdp([rhob;0;pb],p_star);  %R-S-R
    end
end

end