function J = cal_J(W_1,p1,sub)
%Calculate the Riemann invariants
p   = W_1(3,:);

n = 80;
t = linspace(0, 1, n)'; 
P = t * (p1 - p) + p;

rho_star = zeros(n,1);
c_star = zeros(n,1);

for i = 1:n
    rho_star(i) = cal_rho_isentropic(W_1,P(i),sub);
    [c_star(i)]= cal_eq_parameter(P(i),rho_star(i),sub);
end

y = 1./(rho_star.*c_star);
dp = P(2) - P(1);
J = -(dp/3) .* (y(1,:) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end));

end
