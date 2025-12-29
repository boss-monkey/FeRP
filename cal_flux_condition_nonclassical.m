%This program calculates the flow variables on each interface.
function flux_condition = cal_flux_condition_nonclassical(W_L,W_R,star,type_L,States_L)


p_star = star(1,:);
L = star(2,:); R = star(3,:);


rho_L = W_L(1,:); 
u_L = W_L(2,:); u_R = W_R(2,:);
p_L= W_L(3,:);
%%
stateL_a = States_L.a;
stateL_b = States_L.b;
stateL_c = States_L.c; 
stateL_d = States_L.d;

pa_L = stateL_a(1);
rhoa_L = stateL_a(2);
ca_L = stateL_a(3);
pb_L = stateL_b(1);
rhob_L = stateL_b(2);
pc_L = stateL_c(1);
pd_L = stateL_d(1);
W_b_L = W_L;
        
     if type_L == 1
        r_L = 1;
    else
        if p_star > pd_L
            r_L = 1;
        elseif p_star <= pd_L &&  p_star > pb_L
            r_L = 2;
        else
            r_L = 3;
        end
     end

        rho_star_R = cal_rho_hugoniot(W_R, p_star);
        u_star = ((u_L+L)+(u_R-R))/2;      

        c_L = cal_eq_parameter(p_L,rho_L);

        if r_L == 1
            rho_star_L = cal_rho_isentropic(W_L,p_star);
            c_star_L = cal_eq_parameter(p_star,rho_star_L);
            S_L_tail = u_star - c_star_L;
            S_RS_L = S_L_tail;
        elseif r_L == 2
            [p_r_L, rho_r_L, j_r_L] = cal_RS_wave(W_L, p_star, pc_L);
            JL1 = cal_J(W_L,p_r_L);
            u_r_L = u_L + JL1;
            W_r_L = [rho_r_L; u_r_L; p_r_L];
            rho_star_L = cal_rho_isentropic(W_r_L,p_star);
            S_RS_L = u_r_L - j_r_L/rho_r_L;
            S_L_tail = S_RS_L;
        else
            JL1 = cal_J(W_L,pa_L);
            JL2 = (pa_L-pb_L)./(rhoa_L*ca_L);
            ua_L = u_L + JL1;
            ub_L = u_L + JL1 + JL2;
            W_b_L = [rhob_L; ub_L; pb_L];
            rho_star_L = cal_rho_isentropic(W_b_L,p_star);
            c_star_L = cal_eq_parameter(p_star,rho_star_L); 
            S_L_tail = ua_L - ca_L;
            S_RS_L = u_star - c_star_L;

        end

        S_L_head = u_L - c_L;
        S_R = cal_shock_speed(W_R, -1, p_star);
        
       W1_L =  W_L;
       W2_L =  W_L;
       
       x_ref1_L =  S_L_tail;
       x_ref2_L =  S_L_tail;     


       N = 1e4;        

        [~,~,s0_L,~,~,~]= cal_eq_parameter(p_L,rho_L);
        [~,pc] = cal_critical_parameter();
        t = linspace(0, 1, N)'; 
        p1_L = t * (p_star - p_L) + p_L;
        rho1_L = zeros(N,1);
        rho1_L(1) = cal_rho_isentropic0(s0_L, rho_L, p1_L(1));
        for i = 2:N
            rho1_L(i) = cal_rho_isentropic0(s0_L, rho1_L(i-1), p1_L(i));
        end
       for i = 1:N
         if p1_L(i) < pc*0.999
             [~,rho_v,rho_l,~, ~] = cal_T_saturation(p1_L(i));
              if rho1_L(i) > rho_v && rho1_L(i)<rho_l
                  break
              end
         end
       end
       if i ~= N
           c_eq1_L = cal_eq_parameter(p1_L(i-1),rho1_L(i-1));
           c_eq2_L = cal_eq_parameter(p1_L(i),rho1_L(i));
           J1_L = cal_J(W_L,p1_L(i-1));
           J2_L = cal_J(W_L,p1_L(i));
           x_ref1_L = u_L - c_eq1_L + J1_L;
           x_ref2_L = u_L - c_eq2_L + J2_L;
           u1_L = c_eq1_L+x_ref1_L;
           u2_L = c_eq2_L+x_ref2_L;

           W1_L =  [rho1_L(i-1); u1_L ;p1_L(i-1)];
           W2_L =  [rho1_L(i); u2_L ;p1_L(i)];
     
       end

       S = [S_L_head,S_L_tail,S_R,S_RS_L]; 
       
       x_L = [x_ref1_L,x_ref2_L];
       W_star_L =  [rho_star_L;u_star;p_star];
       W_star_R =  [rho_star_R;u_star;p_star];

         
    flux_condition = struct();

    flux_condition.S = S;
    flux_condition.x_L = x_L;
    flux_condition.W_star_L = W_star_L;
    flux_condition.W_star_R = W_star_R;
    flux_condition.W1_L = W1_L;
    flux_condition.W2_L = W2_L;
    flux_condition.W_b_L = W_b_L;


       
end
  
  
  
  
  