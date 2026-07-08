function flux = cal_exact_flux(W_L,W_R,x_ref,flux_condition)
%This program calculates the flow variables on each interface.

S = flux_condition.S;
x_L = flux_condition.x_L;
W_star_L = flux_condition.W_star_L;
W_star_R = flux_condition.W_star_R;
W1_L = flux_condition.W1_L;
W2_L = flux_condition.W2_L;
W_b_L = flux_condition.W_b_L;

   
S_L_head = S(1);
S_L_tail = S(2);
S_R = S(3);
S_RS_L =  S(4);

x_ref1_L = x_L(1);
x_ref2_L = x_L(2);

u_star = W_star_L(2);

    if S_L_head>=x_ref                          % Interface is located in the left undisturbed region
             W = W_L;                           
    elseif S_L_head<x_ref && x_ref<S_L_tail     % Interface is located in the left rarefaction wave
       if x_ref <= x_ref1_L
            W = cal_wave_rho(W_L,1, x_ref, 0);
           elseif x_ref1_L <x_ref && x_ref<= x_ref2_L
            W =  W1_L;
              disp('The interface is located within the split rarefaction wave');
           else
            W = cal_wave_rho(W2_L,1, x_ref, 0);   
       end
    elseif S_L_tail<x_ref && x_ref<S_RS_L 
       W = cal_wave_rho(W_b_L,1, x_ref, 0);
    elseif S_RS_L<x_ref &&  x_ref<S_R   
        if u_star >= x_ref                      % Interface is at the left contact discontinuity
            W =W_star_L;
        else                                    % Interface is at the right contact discontinuity
            W = W_star_R;
        end
    elseif x_ref>=S_R
            W = W_R;                           %Interface is located in the right undisturbed region
    end
    flux = W;

end