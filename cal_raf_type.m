function [type,States] = cal_raf_type(W,p_star)

rho = W(1,:);
p   = W(3,:);

if nargin < 2
   p_star = 200;
end

[~,~,s0,~,~,~]= cal_eq_parameter(p,rho);
[~,pc] = cal_critical_parameter();

N = 2e4;
P = linspace(p, p_star, N);
RHO = zeros(1,N);
RHO(1) = cal_rho_isentropic0(s0, rho, p);

P = [P(1),P];
RHO = [RHO(1),RHO];
disp('————Start calculating the isentropic line boundary————');

for i = 3:N+1

    RHO(i) = cal_rho_isentropic0(s0, RHO(i-1), P(i));

 if P(i) < pc*0.999
     [~,rho_vj,rho_lj,~, ~] = cal_T_saturation(P(i-1));
     [~,rho_vi,~,~, ~] = cal_T_saturation(P(i));
      if RHO(i-1) > rho_vj && RHO(i-1)<rho_lj && RHO(i) < rho_vi          
          break
      end
 end
end
    
 if i == N+1
     disp('Isentropic line boundary calculation complete.');
     type = 1;
     state_a = [0,0,0];
     state_b = [0,0,0];
     state_c = [0,0];
     state_d = [0,0];
     
 else
     disp('Isentropic line boundary calculation complete.');
     type = 2;
     disp('————Start calculating dual-sonic shock————');
     state_c = [P(i-1),RHO(i-1)];
     state_d = [P(i),RHO(i)];
     Pa_guess = P(i-1);
     Pb_guess = P(i);
     [state_a, state_b] = cal_RSR_wave(W, Pa_guess, Pb_guess);
     disp(' '); 
     
 end

States = struct();
States.a = state_a;
States.b = state_b;
States.c = state_c; 
States.d = state_d;
    
end