function [type,States] = cal_raf_type(W, sub,p_star)
%Rarefaction wave types pre-computation
p   = W(3,:);

if nargin < 3
   p_star = 1e4;
end

[~,pc] = cal_critical_parameter(sub);

N = 1e3;
P = linspace(p, p_star, N);
RHO = zeros(1,N);
RHO(1) =  cal_rho_isentropic(W,p,sub);

P = [P(1),P];
RHO = [RHO(1),RHO];
disp('————Start calculating the isentropic line boundary————');

for i = 3:N+1

    RHO(i) = cal_rho_isentropic(W,P(i),sub);

 if P(i) < pc*0.999
     [~,rho_vj,rho_lj,~, ~] = cal_T_saturation(P(i-1),sub);
     [~,rho_vi,~,~, ~] = cal_T_saturation(P(i),sub);
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
     [state_a, state_b] = cal_RSR_wave(W, Pa_guess, Pb_guess,sub);
     disp(' '); 
     
 end

States = struct();
States.a = state_a;
States.b = state_b;
States.c = state_c; 
States.d = state_d;
    
end