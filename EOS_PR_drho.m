function dpdrho = EOS_PR_drho(rho,T)

[a,b,M,~,~] = cal_PR(T);
V = M./rho;
Ru = 8.31443;
dpdV = -Ru.*T./(V - b).^2 + 2.*a.*(V + b)./(V.^2 + 2.*b.*V - b.^2).^2;
dpdrho = -M./rho.^2.*dpdV;
end






