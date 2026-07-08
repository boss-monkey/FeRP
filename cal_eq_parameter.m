function [c_mix,e_mix,s_mix,T,alpha,beta]= cal_eq_parameter(p,rho,sub)
%Compute equilibrium fluid properties: speed of sound, internal energy, entropy, and temperature. 
% 'Sub' denotes the substance (default: $n$-dodecane on the left, nitrogen on the right).
     [~,pc,~,~,M] = cal_critical_parameter(sub);
     N = length(p);
     V = M./rho;
     T =  cal_T_EOS_PR(p,rho,sub);
     c_mix = cal_c_rho_T(rho,T,sub);
     e_mix = cal_energy(rho,T,sub);
     s_mix = cal_entropy(rho,T,sub);
     alpha = zeros(size(rho));
     beta = zeros(size(rho)); 
     for i = 1:N
         if p(i) < pc*0.999
             [Tsat,rho_v,rho_l,V_v, V_l] = cal_T_saturation(p(i),sub);
              if rho(i) > rho_v && rho(i)<rho_l
                 alpha(i) = (rho(i)-rho_l)/(rho_v-rho_l);
                
                 beta(i) = (V(i)-V_l)/(V_v-V_l);
                 T(i) = Tsat;
         
                 e_l = cal_energy(rho_l,Tsat,sub);
                 e_v = cal_energy(rho_v,Tsat,sub);
                 e_mix(i) =  beta(i)*e_v + (1- beta(i))*e_l; 
  
                 
                 s_l = cal_entropy(rho_l,Tsat,sub);
                 s_v = cal_entropy(rho_v,Tsat,sub);
                 s_mix(i) =  beta(i)*s_v + (1- beta(i))*s_l; 
                 
                 h_l = e_l + p(i)./rho_l;
                 h_v = e_v + p(i)./rho_v;
                 dTdP = Tsat.*(1./rho_v - 1./rho_l)./(h_v-h_l);
                 dpdT_l = EOS_PR_dT(rho_l,Tsat,sub);
                 dpdT_v = EOS_PR_dT(rho_v,Tsat,sub);
                 dpdrho_l = EOS_PR_drho(rho_l,Tsat,sub);
                 dpdrho_v = EOS_PR_drho(rho_v,Tsat,sub);
                 drholdP = (1- dpdT_l*dTdP)/dpdrho_l;
                 drhovdP = (1- dpdT_v*dTdP)/dpdrho_v;
                 
                 cp_l =  cal_capacity_p(rho_l,Tsat,sub);
                 cp_v =  cal_capacity_p(rho_v,Tsat,sub);

                 [dsdrho_T_l,dsdT_rho_l] = cal_dsdrho_dT(rho_l,Tsat,sub);
                 [dsdrho_T_v,dsdT_rho_v] = cal_dsdrho_dT(rho_v,Tsat,sub);
                 
                 dsdP_l = dsdrho_T_l.*drholdP +dsdT_rho_l.*dTdP;
                 dsdP_v = dsdrho_T_v.*drhovdP +dsdT_rho_v.*dTdP;
                 cs = Tsat.*(alpha(i)*rho_v/cp_v*dsdP_v^2 + (1-alpha(i))*rho_l/cp_l*dsdP_l^2);
                 
                 c_l = cal_c_rho_T(rho_l,Tsat,sub);
                 c_v = cal_c_rho_T(rho_v,Tsat,sub);

                 c_mix(i) = (rho(i)*(alpha(i)/rho_v/c_v^2 + (1-alpha(i))/rho_l/c_l^2+cs))^(-0.5); 
              elseif rho(i) < rho_v
                   alpha(i) = 1;
                   beta(i) = 1;
             end
         end
  
     end
end