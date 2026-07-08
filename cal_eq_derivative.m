function [dcdrho_p,dcdp_rho,dedrho_p,dedp_rho,dsdrho_p,dsdp_rho]= cal_eq_derivative(p,rho,sub)
%Compute equilibrium fluid derivatives (speed of sound, internal energy, entropy, and temperature). 
% 'Sub' denotes the substance; default is $n$-dodecane on the left and nitrogen on the right.
     [~,pc,~,~,~] = cal_critical_parameter(sub);
     N = length(p);
     T =  cal_T_EOS_PR(p,rho,sub);
     [dedrho_p,dedp_rho] = cal_dedrho_dp(rho,T,sub);
     [dsdrho_p,dsdp_rho] = cal_dsdrho_dp(rho,T,sub);
     [dcdrho_p,dcdp_rho] = cal_dcdrho_dp(rho,T,sub);
     for i = 1:N
         if p(i) < pc*0.999
             [Tsat,rho_v,rho_l,~, ~] = cal_T_saturation(p(i),sub);
              if rho(i) > rho_v && rho(i)<rho_l
                 pi = p(i); rhoi = rho(i);
                 alpha = (rhoi-rho_l)/(rho_v-rho_l); 
                 v_v = 1./rho_v; v_l = 1./rho_l; vi = 1./rhoi;
                 beta =  (vi-v_l)/(v_v-v_l);
                 T(i) = Tsat;
         
                 e_l = cal_energy(rho_l,Tsat,sub);
                 e_v = cal_energy(rho_v,Tsat,sub);
  
                 
                 s_l = cal_entropy(rho_l,Tsat,sub);
                 s_v = cal_entropy(rho_v,Tsat,sub);
                 
 %%
                 h_l = e_l + pi/rho_l;
                 h_v = e_v + pi./rho_v;
                 dh = (h_v-h_l);
                 dv = (v_v - v_l);
                 dTdP = Tsat.*dv./dh;

                 drhodP_l = cal_drhodP_sat(rho_l, Tsat, dTdP,sub);
                 drhodP_v = cal_drhodP_sat(rho_v, Tsat, dTdP,sub);
                 
                 c_l = cal_c_rho_T(rho_l,Tsat,sub);
                 c_v = cal_c_rho_T(rho_v,Tsat,sub);
                 cp_l =  cal_capacity_p(rho_l,Tsat,sub);
                 cp_v =  cal_capacity_p(rho_v,Tsat,sub);
                 C_v = 1/rho_v/c_v^2; 
                 C_l = 1/rho_l/c_l^2;
                 Kw = rhoi.*(alpha*(C_v-C_l) + C_l);
                 dKwdrho_p = Kw./rhoi+rhoi*(C_v-C_l)/(rho_v-rho_l);
                 cw = (Kw)^(-0.5); 
                 
                 [dcdrho_l,dcdT_l] = cal_dcdrho_dT(rho_l,Tsat,sub);
                 [dcdrho_v,dcdT_v] = cal_dcdrho_dT(rho_v,Tsat,sub);
                 
                 dcl2dP = 2*c_l*(dcdrho_l*drhodP_l + dcdT_l*dTdP);
                 dcv2dP = 2*c_v*(dcdrho_v*drhodP_v + dcdT_v*dTdP);
                 
                 dCldP = -C_l^2*(drhodP_l*c_l^2+rho_l*dcl2dP);
                 dCvdP = -C_v^2*(drhodP_v*c_v^2+rho_v*dcv2dP);
                 
                 dalphadp = (rho_l * drhodP_v - rho_v*drhodP_l + rhoi*(drhodP_l-drhodP_v))/(rho_l-rho_v)^2;     
                 dKwdp_rho = rhoi*(dalphadp * (C_v-C_l) + alpha*(dCvdP - dCldP) + dCldP);

                
                
%%
                 dbetadp = rho_v./rhoi*dalphadp + alpha./rhoi.*drhodP_v;

                 demdrho_p = -(e_v-e_l)./((v_v-v_l).*rhoi.^2);
                 dsmdrho_p = -(s_v-s_l)./((v_v-v_l).*rhoi.^2);


                 cv_l = cal_capacity(rho_l,Tsat,sub);
                 cv_v =  cal_capacity(rho_v,Tsat,sub);
          
                 [dsdrho_T_l,dsdT_rho_l] = cal_dsdrho_dT(rho_l,Tsat,sub);
                 [dsdrho_T_v,dsdT_rho_v] = cal_dsdrho_dT(rho_v,Tsat,sub);
                 
                 dsdP_l = dsdrho_T_l.*drhodP_l +dsdT_rho_l.*dTdP;
                 dsdP_v = dsdrho_T_v.*drhodP_v +dsdT_rho_v.*dTdP;
                 
                 dedT_rho_l = cv_l;
                 dedT_rho_v = cv_v;
                 
                 dedrho_T_l = cal_dedrho_T(rho_l,Tsat,sub);
                 dedrho_T_v = cal_dedrho_T(rho_v,Tsat,sub);
                 
                 dedP_l = dedrho_T_l.*drhodP_l +dedT_rho_l.*dTdP;
                 dedP_v = dedrho_T_v.*drhodP_v +dedT_rho_v.*dTdP;
                 
                 demdp_rho = beta.*dedP_v + (1-beta).*dedP_l+(e_v-e_l).*dbetadp;
                 dsmdp_rho = beta.*dsdP_v + (1-beta).*dsdP_l+(s_v-s_l).*dbetadp;

                 dedrho_p(i) = demdrho_p;
                 dedp_rho(i) = demdp_rho;
                 
                 dsdrho_p(i) = dsmdrho_p;
                 dsdp_rho(i) = dsmdp_rho;

                 %%
                 Kp = rhoi*Tsat.*(alpha*rho_v/cp_v*dsdP_v^2 + (1-alpha)*rho_l/cp_l*dsdP_l^2);
                 
                 f = Tsat*rho_v/cp_v*dsdP_v^2;
                 g = Tsat*rho_l/cp_l*dsdP_l^2;
                 
                 dvdP_l = -1/rho_l^2* drhodP_l;
                 dvdP_v = -1/rho_v^2* drhodP_v;
                 dhdP_l = dedP_l +v_l + pi*dvdP_l;
                 dhdP_v = dedP_v +v_v + pi*dvdP_v;
                 
                 d2TdP2 = 1/dh^2 *((dTdP*dv + Tsat*(dvdP_v - dvdP_l))*dh - Tsat*dv*(dhdP_v -dhdP_l));
                 
                 d2sdP2_l = cal_d2sdP2_sat(rho_l,Tsat,dTdP,d2TdP2,sub);
                 d2sdP2_v = cal_d2sdP2_sat(rho_v,Tsat,dTdP,d2TdP2,sub);
                 
                 
                 [dcpdrho_l, dCpdT_l] = cal_dCpdrho_dT(rho_l, Tsat, sub);
                 [dcpdrho_v, dCpdT_v] = cal_dCpdrho_dT(rho_v, Tsat, sub);
                 dcpdP_l = dcpdrho_l.*drhodP_l +dCpdT_l*dTdP;
                 dcpdP_v = dcpdrho_v.*drhodP_v +dCpdT_v*dTdP;
                 
                 dfdP = f*(dTdP/Tsat + drhodP_v/rho_v - dcpdP_v/cp_v +2*d2sdP2_v/dsdP_v);
                 dgdP = g*(dTdP/Tsat + drhodP_l/rho_l - dcpdP_l/cp_l +2*d2sdP2_l/dsdP_l);
                 
                 dKpdp_rho = rhoi.*((f-g)*dalphadp +alpha*(dfdP - dgdP) +dgdP);
                 
                 c_eq = (Kw + Kp)^(-0.5); 
                 dKpdrho_p = Kp/rho(i) +rhoi*(f-g)/(rho_v-rho_l);
                 
                 dcdrho_p(i) = -c_eq^3 * (dKwdrho_p + dKpdrho_p)/2;
                 dcdp_rho(i) = -c_eq^3 * (dKwdp_rho + dKpdp_rho)/2;

             end
         end
  
     end
end