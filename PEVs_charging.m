% Numerical example for the manuscript
% "Perturbed Distributed Online Optimization", by
% Yanxu Su, Qingling Wang and Changyin Sun, 2021

clear all
close all
clc
load Init_data.mat
%%
rng(238);
%Numerical example problem definition
% Dimension of x
n = 24;

% Dimension of constraints
m = 48;

% Time Horizon
T = 50;
epson = 0.5;
epson1 = 0.3;
epson2 = 0.5;
epson3 = 0.7;
epson4 = 0.9;
epsilon = 0.01;

sigma1 = 2;
sigma2 = 1;

% Time-interval between strongly connected digrapghs
Q = 4;

% Time-varying networks description
cnt =0;
N = 50; %Number of subsystems

% Initialization of variables
x = zeros(n,T,N);
tilde_x = zeros(n,T,N);
tilde_x_c = zeros(n,T,N);
mu = zeros(m,T,N);
tilde_mu = zeros(m,T,N);
y = zeros(n,T,N);
x(:,1,:) = 0.5*rand(n,N);
x0(:,1,:) = 5*rand(n,N);

Ci = 10*rand(n,T,N);

Pi_min = 0.1*ones(n,N);
Pi_max = 1.5*ones(n,N);
Ei_min = ones(N,1);
Ei_max = 8+(16-8)*rand(N,1);
Ei_ini = (0.2+(0.5-0.2)*rand(N,1)).*Ei_max;
Ei_ref = (0.55+(0.8-0.55)*rand(N,1)).*Ei_max;

x_c_e_1 = zeros(n,T,N);
tilde_x_c_e_1 = zeros(n,T,N);
mu_c_e_1 = zeros(m,T,N);
tilde_mu_c_e_1 = zeros(m,T,N);
x_c_e_1(:,1,:) = x(:,1,:);
x_c_e_2 = zeros(n,T,N);
tilde_x_c_e_2 = zeros(n,T,N);
mu_c_e_2 = zeros(m,T,N);
tilde_mu_c_e_2 = zeros(m,T,N);
x_c_e_2(:,1,:) = x(:,1,:);
x_c_e_3 = zeros(n,T,N);
tilde_x_c_e_3 = zeros(n,T,N);
mu_c_e_3 = zeros(m,T,N);
tilde_mu_c_e_3 = zeros(m,T,N);
x_c_e_3(:,1,:) = x(:,1,:);
x_c_e_4 = zeros(n,T,N);
tilde_x_c_e_4 = zeros(n,T,N);
mu_c_e_4 = zeros(m,T,N);
tilde_mu_c_e_4 = zeros(m,T,N);
x_c_e_4(:,1,:) = x(:,1,:);

x_c = zeros(n,T,N);
x_c(:,1,:) = x(:,1,:);
q_c = zeros(m,T,N);
tilde_q = zeros(m,T,N);
c = 0.5;
kappa = 0.5;
zeta0 = 0.015+(0.075-0.015)*rand(N,1);
zeta = (1-zeta0)/3;

% Model description
for i=1:N
    H_i(:,:,i) = [eye(n); -eye(n)];% -5 + 10* rand(m,n);
end

zeta_1 = 1;
zeta_2 = 5;

V = zeros(T,1);
V_c_e_1 = zeros(T,1);
V_c_e_2 = zeros(T,1);
V_c_e_3 = zeros(T,1);
V_c_e_4 = zeros(T,1);
V_c = zeros(T,1);

A = geneDbStoch(n);

d = N*ones(m,1);
l = rand(1);
delta(:,1) = N*(0/T)*rand(m,1);
for j=1:m
    delta_T_u(j,1) = max(delta(j,1));
end

switch N
    case 50
        S = Q_connected(1:3,:);
    case 100
        S = Q_connected(4:6,:);
    case 200
        S = Q_connected(7:9,:);
end
switch Q
    case 4
        S = Q_connected(1,:);
    case 8
        S = Q_connected(2,:);
    case 16
        S = Q_connected(3,:);
end

q = S(3);
r = S(4);
Q_cnt = 0;

for iter = 1:T
    W(:,:,iter) = eye(N);
    Digrapgh = genRandomStronglyConnectedGraph(q);
    if q+Q_cnt*r < N
        for i=1+Q_cnt*r:q+Q_cnt*r
            for j=1+Q_cnt*r:q+Q_cnt*r
                if (j ~= i)&&(Digrapgh(i-Q_cnt*r,j-Q_cnt*r)>0)
                    W(i,j,iter) = 1/q;
                elseif j == i
                    for k = 1+Q_cnt*r:q+Q_cnt*r
                        if (Digrapgh(i-Q_cnt*r,k-Q_cnt*r)>0)&&(k~=i)
                            cnt = cnt + 1;
                        end
                    end
                    W(i,j,iter) = 1-cnt/q;
                    cnt = 0;
                else
                    W(i,j,iter) = 0;
                end
            end
        end
    else
        for i=1+Q_cnt*r:N
            for j=1+Q_cnt*r:N
                if (j ~= i)&&(Digrapgh(i-Q_cnt*r,j-Q_cnt*r)>0)
                    W(i,j,iter) = 1/q;
                elseif j == i
                    for k = 1+Q_cnt*r:N
                        if (Digrapgh(i-Q_cnt*r,k-Q_cnt*r)>0)&&(k~=i)
                            cnt = cnt + 1;
                        end
                    end
                    W(i,j,iter) = 1-cnt/q;
                    cnt = 0;
                else
                    W(i,j,iter) = 0;
                end
            end
        end
        for i=1:q+Q_cnt*r-N
            for j=1:q+Q_cnt*r-N
                if (j ~= i)&&(Digrapgh(i,j)>0)
                    W(i,j,iter) = 1/q;
                elseif j == i
                    for k = 1:q+Q_cnt*r-N
                        if (Digrapgh(i,k)>0)&&(k~=i)
                            cnt = cnt + 1;
                        end
                    end
                    W(i,j,iter) = 1-cnt/q;
                    cnt = 0;
                else
                    W(i,j,iter) = 0;
                end
            end
        end
    end
    if Q_cnt < Q-1
        Q_cnt = Q_cnt +1;
    else
        Q_cnt = 0;
    end
    
    for i = 1:N
        for j = 1:N
            if (W(i,j,iter) > 0)
                tilde_mu(:,iter,i) = tilde_mu(:,iter,i)+ W(i,j,iter)*mu(:,iter,j);
            end
        end
        rho = 1/(iter^epson);
        Df = Ci(:,iter,i);
        %g(:,iter,i) = H_i(:,:,i)*x(:,iter,i)-d(i,:)'+delta(:,iter,i);
        x(:,iter+1,i) = myArgMin(Df, tilde_mu(:,iter,i), H_i(:,:,i), d/N,...
            rho, x(:,iter,i), n, zeta(i), Ei_ini(i), Ei_ref(i), Ei_min(i), Ei_max(i));
        g(:,iter+1,i) = H_i(:,:,i)*x(:,iter+1,i)-d/N;
    end
    l(iter+1) = rand(1);
    delta(:,iter+1) = ((iter+1)*0.4/T)*N*rand(m,1);
    for i = 1:N
        mu(:,iter+1,i) = max(tilde_mu(:,iter,i)+rho*(g(:,iter+1,i)+delta(:,iter+1)/N)/(2*sigma2),0);
        %myArgMax(tilde_mu(:,iter,i), g(:,iter+1,i), delta(:,iter+1)/N, rho, m);
        f(i,iter) = Ci(:,iter,i)'*x(:,iter,i);
        
        for j=1:m
            delta_T_u(j,iter) = max(delta(j,:));
        end
        V(iter,1) = V(iter,1) + norm(max(H_i(:,:,i)*x(:,iter,i)-d/N+delta(:,iter)/N,0))/iter;
    end
    Reg1(iter) = sum(sum(f))/iter;
    
    % Compared with different epson
    %     for i = 1:N
    %         %z(:,iter+1,i) = sum(x(:,iter+1,i));
    %         %tilde_z(:,iter+1,i) = sum(y(:,iter+1,i));
    %         for j = 1:N
    %             if (W(i,j,iter) > 0)
    %                 %tilde_x_c_e_1(:,iter,i) = tilde_x_c_e_1(:,iter,i) + W(i,j)*x_c_e_1(:,iter,j);
    %                 tilde_mu_c_e_1(:,iter,i) = tilde_mu_c_e_1(:,iter,i)+ W(i,j,iter)*mu_c_e_1(:,iter,j);
    %             end
    %         end
    %         rho1 = 1/(iter^epson1);
    %         Df = Ci(:,iter,i);%zeta_1*pi(:,iter,i)+2*zeta_2*(x_c_e_1(:,iter,i)-y(:,iter,i));
    %         %g(:,iter,i) = H_i(:,:,i)*x(:,iter,i)-d(i,:)'+delta(:,iter,i);
    %         x_c_e_1(:,iter+1,i) = myArgMin(Df, tilde_mu_c_e_1(:,iter,i), H_i(:,:,i), d/N,...
    %             rho1, x_c_e_1(:,iter,i), n, zeta(i), Ei_ini(i), Ei_ref(i), Ei_min(i), Ei_max(i));
    %         g_c_e_1(:,iter+1,i) = H_i(:,:,i)*x_c_e_1(:,iter+1,i)-d/N;
    %
    %         mu_c_e_1(:,iter+1,i) = myArgMax(tilde_mu_c_e_1(:,iter,i), ...
    %             g_c_e_1(:,iter+1,i), delta(:,iter+1)/N, rho1, m);
    %         f_c_e_1(i,iter) = Ci(:,iter,i)'*x_c_e_1(:,iter,i);
    %         %f_star(i,iter) = zeta_1*pi(:,iter,i)'*x_T(:,i)+zeta_2*(x_T(:,i)-y(:,iter,i))'*(x_T(:,i)-y(:,iter,i));
    %         V_c_e_1(iter,1) = V_c_e_1(iter,1) + norm(max(H_i(:,:,i)*x_c_e_1(:,iter,i)-d/N+delta(:,iter)/N,0))/iter;
    %     end
    %     Reg_c_e_1(iter) = sum(sum(f_c_e_1))/iter;
    %
    %     for i = 1:N
    %         for j = 1:N
    %             if (W(i,j,iter) > 0)
    %                 %tilde_x_c_e_2(:,iter,i) = tilde_x_c_e_2(:,iter,i) + W(i,j)*x_c_e_2(:,iter,j);
    %                 tilde_mu_c_e_2(:,iter,i) = tilde_mu_c_e_2(:,iter,i)+ W(i,j,iter)*mu_c_e_2(:,iter,j);
    %             end
    %         end
    %         rho2 = 1/(iter^epson2);
    %         Df = Ci(:,iter,i);
    %         %g(:,iter,i) = H_i(:,:,i)*x(:,iter,i)-d(i,:)'+delta(:,iter,i);
    %         x_c_e_2(:,iter+1,i) = myArgMin(Df, tilde_mu_c_e_2(:,iter,i), H_i(:,:,i), d/N,...
    %             rho2, x_c_e_2(:,iter,i), n, zeta(i), Ei_ini(i), Ei_ref(i), Ei_min(i), Ei_max(i));
    %         g_c_e_2(:,iter+1,i) = H_i(:,:,i)*x_c_e_2(:,iter+1,i)-d/N;
    %         mu_c_e_2(:,iter+1,i) = myArgMax(tilde_mu_c_e_2(:,iter,i), ...
    %             g_c_e_2(:,iter+1,i), delta(:,iter+1)/N, rho2, m);
    %         f_c_e_2(i,iter) = Ci(:,iter,i)'*x_c_e_2(:,iter,i);
    %
    %         V_c_e_2(iter,1) = V_c_e_2(iter,1) + norm(max(H_i(:,:,i)*x_c_e_2(:,iter,i)-d/N+delta(:,iter)/N,0))/iter;
    %     end
    %     Reg_c_e_2(iter) = sum(sum(f_c_e_2))/iter;
    %
    %     for i = 1:N
    %         for j = 1:N
    %             if (W(i,j,iter) > 0)
    %                 %tilde_x_c_e_3(:,iter,i) = tilde_x_c_e_3(:,iter,i) + W(i,j)*x_c_e_3(:,iter,j);
    %                 tilde_mu_c_e_3(:,iter,i) = tilde_mu_c_e_3(:,iter,i)+ W(i,j,iter)*mu_c_e_3(:,iter,j);
    %             end
    %         end
    %         rho3 = 1/(iter^epson3);
    %         Df = Ci(:,iter,i);
    %         %g(:,iter,i) = H_i(:,:,i)*x(:,iter,i)-d(i,:)'+delta(:,iter,i);
    %         x_c_e_3(:,iter+1,i) = myArgMin(Df, tilde_mu_c_e_3(:,iter,i), H_i(:,:,i), d/N,...
    %             rho3, x_c_e_3(:,iter,i), n, zeta(i), Ei_ini(i), Ei_ref(i), Ei_min(i), Ei_max(i));
    %         g_c_e_3(:,iter+1,i) = H_i(:,:,i)*x_c_e_3(:,iter+1,i)-d/N;
    %         mu_c_e_3(:,iter+1,i) = myArgMax(tilde_mu_c_e_3(:,iter,i), ...
    %             g_c_e_3(:,iter+1,i), delta(:,iter+1)/N, rho3, m);
    %         f_c_e_3(i,iter) = Ci(:,iter,i)'*x_c_e_3(:,iter,i);
    %
    %         V_c_e_3(iter,1) = V_c_e_3(iter,1) + norm(max(H_i(:,:,i)*x_c_e_3(:,iter,i)-d/N+delta(:,iter)/N,0))/iter;
    %     end
    %     Reg_c_e_3(iter) = sum(sum(f_c_e_3))/iter;
    %
    %     for i = 1:N
    %         for j = 1:N
    %             if (W(i,j,iter) > 0)
    %                 %tilde_x_c_e_3(:,iter,i) = tilde_x_c_e_3(:,iter,i) + W(i,j)*x_c_e_3(:,iter,j);
    %                 tilde_mu_c_e_4(:,iter,i) = tilde_mu_c_e_4(:,iter,i)+ W(i,j,iter)*mu_c_e_4(:,iter,j);
    %             end
    %         end
    %         rho4 = 1/(iter^epson4);
    %         Df = Ci(:,iter,i);
    %         %g(:,iter,i) = H_i(:,:,i)*x(:,iter,i)-d(i,:)'+delta(:,iter,i);
    %         x_c_e_4(:,iter+1,i) = myArgMin(Df, tilde_mu_c_e_4(:,iter,i), H_i(:,:,i), d/N,...
    %             rho4, x_c_e_4(:,iter,i), n, zeta(i), Ei_ini(i), Ei_ref(i), Ei_min(i), Ei_max(i));
    %         g_c_e_4(:,iter+1,i) = H_i(:,:,i)*x_c_e_4(:,iter+1,i)-d/N;
    %         mu_c_e_4(:,iter+1,i) = myArgMax(tilde_mu_c_e_4(:,iter,i), ...
    %             g_c_e_4(:,iter+1,i), delta(:,iter+1)/N, rho4, m);
    %         f_c_e_4(i,iter) = Ci(:,iter,i)'*x_c_e_4(:,iter,i);
    %
    %         V_c_e_4(iter,1) = V_c_e_4(iter,1) + norm(max(H_i(:,:,i)*x_c_e_4(:,iter,i)-d/N+delta(:,iter)/N,0))/iter;
    %     end
    %     Reg_c_e_4(iter) = sum(sum(f_c_e_4))/iter;
    for i =1:N
        delta_T_l(:,iter) = sum(delta,2)/iter;
    end
    
end

%% Best Minimizers (X_T)

cvx_begin
variable x_B(n,N)
variable e(n,N)
L = 0;

G1 = zeros(m,1);
G2 = zeros(m,1);
w  = zeros(m,T);

for i = 1:N    
    L = L+ sum(x_B(:,i)'*Ci(:,:,i));
    G1 = G1 + H_i(:,:,i)*x_B(:,i);
end
for t = 1:T
    w(:,t) = (1/(t^epson))*delta_T_l(:,T)+ (1-1/(t^epson))*delta_T_u(:,T);
end
minimize L
subject to
Pi_min <= x_B <= Pi_max;
e(1,:) == Ei_ini';
for jj=1:n-1
    e(jj+1,:) == e(jj,:)+ zeta'.*x_B(jj,:);
    Ei_min' <= e(jj+1,:) <= Ei_max';
end
e(n,:) == Ei_ref';
for t = 1:T
    G1-d+w(:,t) <=0;
end
cvx_end
x_T = x_B;

%% Best Minimizers (X_T_max)
cvx_begin
variable x_B_max(n,N)
variable e_max(n,N)
L_max = 0;

G1_max = zeros(m,1);

for i = 1:N
    L_max = L_max+ sum(x_B_max(:,i)'*Ci(:,:,i));
    G1_max = G1_max + H_i(:,:,i)*x_B_max(:,i);
end

minimize L_max
subject to
Pi_min <= x_B_max <= Pi_max;
e_max(1,:) == Ei_ini';
for jj=1:n-1
    e_max(jj+1,:) == e_max(jj,:)+ zeta'.*x_B_max(jj,:);
    Ei_min' <= e_max(jj+1,:) <= Ei_max';
end
e_max(n,:) == Ei_ref';
G1_max-d+delta_T_l(:,T) <=0;
cvx_end
x_T_max = x_B_max;

%% Best Minimizers (X_T_min)
cvx_begin
variable x_B_min(n,N)
variable e_min(n,N)
L_min = 0;

G1_min = zeros(m,1);

for i = 1:N
    L_min = L_min+ sum(x_B_min(:,i)'*Ci(:,:,i));
    G1_min = G1_min + H_i(:,:,i)*x_B_min(:,i);
end
minimize L_min
subject to
Pi_min <= x_B_min <= Pi_max;
e_min(1,:) == Ei_ini';
for jj=1:n-1
    e_min(jj+1,:) == e_min(jj,:)+ zeta'.*x_B_min(jj,:);
    Ei_min' <= e_min(jj+1,:) <= Ei_max';
end
e_min(n,:) == Ei_ref';
G1_min-d+delta_T_u(:,T) <=0;
cvx_end
x_T_min = x_B_min;


f_star = zeros(N,T);
f_star_max = zeros(N,T);
f_star_min = zeros(N,T);
Reg = zeros(T,1);
Reg_max = zeros(T,1);
Reg_min = zeros(T,1);
for t=1:T
    for i=1:N
        f_star(i,t) = Ci(:,t,i)'*x_T(:,i);
        f_star_max(i,t) = Ci(:,t,i)'*x_T_max(:,i);
        f_star_min(i,t) = Ci(:,t,i)'*x_T_min(:,i);
    end
    Reg(t) = Reg1(t)-sum(sum(f_star))/t;
    Reg_max(t) = Reg1(t)-sum(sum(f_star_max))/t;
    Reg_min(t) = Reg1(t)-sum(sum(f_star_min))/t;
%         Reg_c_1(t) = Reg_c_e_1(t)-sum(sum(f_star))/t;
%         Reg_c_2(t) = Reg_c_e_2(t)-sum(sum(f_star))/t;
%         Reg_c_3(t) = Reg_c_e_3(t)-sum(sum(f_star))/t;
%         Reg_c_4(t) = Reg_c_e_4(t)-sum(sum(f_star))/t;
    %    Reg_c(t) = Reg2(t)-sum(sum(f_star))/t;
end
% for t=1:T
% Reg33(t) = Reg1(t)-Reg2(t);
% end

set(0, 'DefaultAxesFontSize',15)
set(0, 'DefaultTextFontSize',15)
figure(1);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
xlabel('$T$','Interpreter','Latex','FontSize',15);
ylabel('$V(T)/T$','Interpreter','Latex','FontSize',15);
grid on;
hold on;
plot(V,'b','LineWidth',2);
plot(V8,'r--','LineWidth',2);
plot(V16,'g-.','LineWidth',2);
% plot(V_c_e_1,'r--','LineWidth',2);
% plot(V_c_e_2,'g-.','LineWidth',2);
% plot(V_c_e_3,'c:','LineWidth',2);
%plot(V_c,'r--','LineWidth',2);
%plot(V_c2,'r','LineWidth',2);
figure(2);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
xlabel('$T$','Interpreter','Latex','FontSize',15);
ylabel('$R(T)/T$','Interpreter','Latex','FontSize',15);
grid on;
hold on;
plot(Reg,'b','LineWidth',2);
plot(Reg_max,'r--','LineWidth',2);
plot(Reg_min,'g-.','LineWidth',2);

% plot(Reg_c_1,'r--','LineWidth',2);
% plot(Reg_c_2,'g-.','LineWidth',2);
% plot(Reg_c_3,'c:','LineWidth',2);
%plot(Reg_c,'r--','LineWidth',2);
%plot(Reg33,'c:','LineWidth',2);

