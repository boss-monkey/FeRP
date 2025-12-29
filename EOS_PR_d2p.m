function [d2pdrho2,d2pdTdrho,d2pdT2] = EOS_PR_d2p(rho,T)

[a,b,M,dadT,d2adT2] = cal_PR(T);
Ru =8.31443;

V = M./rho;

dpdV = -Ru.*T./(V - b).^2 + 2.*a.*(V + b)./(V.^2 + 2.*b.*V - b.^2).^2;
d2pdV2 = 2.*Ru.*T./(V - b).^3 - 2.*a.*(3.*V.^2 + 6.*b.*V + 5.*b.^2)./(V.^2 + 2.*b.*V - b.^2).^3;
d2pdrho2 = (M^2 ./ rho.^4) .* d2pdV2 + (2.*M ./ rho.^3) .* dpdV;
d2pdTdrho = M./rho.^2.*( Ru./(V-b).^2 - 2*dadT.*(V+b)./(V.^2+2*V.*b-b.^2).^2 ); 
d2pdT2 = -d2adT2./(V.^2+2*V.*b-b.^2);

end