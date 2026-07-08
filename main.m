clearvars
clc
p_0 = 2.0e+06;  p_1 = 1e+05;
rho_0 = 400;  rho_1 = 1;
u_0   =200; u_1   = 0;


[~,~,~,T_0]= cal_eq_parameter(p_0,rho_0, 0);
[~,~,~,T_1]= cal_eq_parameter(p_1,rho_1, 1);
%% Set calculation domain and calculation time
time = 600e-6;                    %time
x1 = 0; x2 = 1; x3 = 0.5;       %range
N = 500;                        %grid number
x = linspace(x1,x2,N);

%% Calculat the Intermediate pressure
W_L = [rho_0; u_0 ;p_0];
W_R = [rho_1; u_1 ;p_1];   %The left and right state

%%
[type_L,States_L] = cal_raf_type(W_L, 0);
if type_L == 1
     disp('The left isentropic line does not pass positively through the gas phase saturation line.');
else
     disp('The left isentropic line passes positively through the gas phase saturation line, successfully calculating the dual-sonic shock.');
end


%%
tic
star = cal_p_star_nonclassical(W_R,W_L,type_L,States_L);
toc
disp('————Start calculating interface flux————'); 
%% Calculate the flow variables on each interface
tic

W_1 = zeros(3,N);
flux_condition = cal_flux_condition_nonclassical(W_L,W_R,star,type_L,States_L);
x_ref =(x-x3)/time;
for i = 1:N
    flux = cal_exact_flux(W_L,W_R,x_ref(i),flux_condition);
    W_1(:,i) = flux;
    disp(['Current grid interface' num2str(i)]); 
end
toc
%% Plot
rho1 = W_1(1,:); u1 = W_1(2,:);
p1   = W_1(3,:);
u_star = flux_condition.W_star_L(2);
[c1,~,~,T1,alpha,beta]= cal_eq_parameter(p1,rho1,0);

for i = 1:N
    if u_star < x_ref(i)  
    [c1(i),~,~,T1(i),alpha(i),beta(i)]= cal_eq_parameter(p1(i),rho1(i),1);
    end
end

result = [x',rho1',u1',p1',T1',c1',alpha'];

set(groot,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex',...
    'defaultTextInterpreter','latex',...
    'DefaultTextFontSize', 16, ...
    'DefaultAxesFontSize', 18, ...
    'DefaultLegendFontSize', 18);

figure(1);
subplot(2,3,1)      % pressure
plot(x,p1,'-r','linewidth',2)
xlabel('x');
ylabel('$p$');
% set(gca, 'YScale', 'log')
grid on
subplot(2,3,2)      % u
plot(x,u1,'-b','linewidth',2)
xlabel('x');
ylabel('$u$');
grid on
subplot(2,3,3)      % temperature
plot(x,T1,'-k','linewidth',2)
xlabel('x');
ylabel('$T$');
grid on
subplot(2,3,4)     % rho
plot(x,rho1,'-m','linewidth',2)
xlabel('x');
ylabel('$\rho$');
grid on
subplot(2,3,5)      % c
plot(x,c1,'-b','linewidth',2)
xlabel('x');
ylabel('$c_{eq}$');
grid on
sgtitle('Riemann problem of Doedecane-N2', 'FontSize', 16);
subplot(2,3,6)      % beta
plot(x,alpha,'-r','linewidth',2)
xlabel('x');
ylabel('$\beta$');
grid on