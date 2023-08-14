clear all, close all;
%table = readmatrix('mixing_system.xlsx');
data2 = load('data2.mat');
t = data2.t;
table = data2.x;
c1 = movmean(table(1,:),10);
c2 = movmean(table(2,:),10);
c3 = movmean(table(3,:),10);

figure
hold on
plot(c1)
plot(c2)
plot(c3)
title('Measured values with noise reduction')
legend('c1','c2','c3')

%h = 1.0;
c01 = 1;
c02 = 2;
V = 50;

N = length(t);
alpha = 1;          % haven't used
lambda = 0.95;      % haven't used

F1_est = zeros(N,1);
F2_est = zeros(N,1);
k_est = zeros(N,1);
F_est = zeros(N,1);
P = zeros(3,1);     % haven't used
%phi = zeros(4,3);
Y = zeros(1,3);

for i = 1:N
    %phi(2,:) = phi(1,:);
    %phi(1,:) = Y;
    %phi(4,:) = phi(3,:);
    %phi(3,:) = [c1(i) c2(i) c3(i)];
    phi = [c1(i) c2(i) c3(i)];  % haven't used

    X = [c01/V 0 -c1(i)*c2(i) -c1(i)/V;...
         0 c02/V -c1(i)*c2(1) -c2(i)/V;...
         0     0 -c1(i)*c2(i) -c3(i)/V];
    if i == 1
        Y = [c1(1)/t(1) c2(1)/t(1) c3(1)/t(1)];
    else
        Y = [(c1(i)-c1(i-1))/(t(i)-t(i-1))...
             (c2(i)-c2(i-1))/(t(i)-t(i-1))...
             (c3(i)-c3(i-1))/(t(i)-t(i-1))];
    end
    
    %{
    P(1) = [1/lambda*(P(1)-(P(1)*phi*transpose(phi)*P(1))/(lambda/alpha)+transpose(phi)*P(1)*phi)];
    P(2) = [1/lambda*(P(2)-(P(2)*phi*transpose(phi)*P(2))/(lambda/alpha)+transpose(phi)*P(2)*phi)];
    P(3) = [1/lambda*(P(3)-(P(3)*phi*transpose(phi)*P(3))/(lambda/alpha)+transpose(phi)*P(3)*phi)];
    %}
    theta = inv(transpose(X)*X)*transpose(X)*transpose(Y);
    F1_est(i) = theta(1);
    F2_est(i) = theta(2);
    k_est(i) = theta(3);
    F_est(i) = theta(4);
end

%{
for i = 2:N
    syms F1 F2 F k
    [F1_ans, F2_ans, F_ans, k_ans] = solve([F1*c01/V-k*c1(i-1)*c2(i-1)-F*c1(i-1)/V-(c1(i)-c1(i-1)),...
        F2*c02/V-k*c1(i-1)*c2(i-1)-F*c2(i-1)/V-(c2(i)-c2(i-1)),  k*c1(i-1)*c2(i-1)-F*c3(i-1)/V-(c3(i)-c3(i-1))],...
        F1+F2-F,[F1 F2 F k]);
    F1_est(i) = F1_ans;
    F2_est(i) = F2_ans;
    F_est(i) = F_ans;
    k_est(i) = k_ans;
end
%}

figure
hold on
plot(F1_est)
plot(F2_est)
plot(F_est)
legend('F1 estimation','F2 estimation','F estimation')
figure
plot(k_est)
legend('k estimation')

%% prova
c1_est = zeros(size(c1));
c2_est = zeros(size(c1));
c3_est = zeros(size(c1));

c1_est(1) = c1(1);
c2_est(1) = c2(1);
c3_est(1) = c3(1);

dt = 0.1;
for i = 2:500
    c1_est(i) = c1_est(i-1) + F1_est(i)*c01/V - k_est(i)*c1_est(i)*c2_est(i) - F_est(i)*c1_est(i)/V;
    c2_est(i) = c2_est(i-1) + F2_est(i)*c02/V - k_est(i)*c1_est(i)*c2_est(i) - F_est(i)*c2_est(i)/V;
    c3_est(i) = k_est(i)*c1_est(i)*c2_est(i) - F_est(i)*c3_est(i)/V;
end

figure
hold on
plot(c1_est)
plot(c2_est)
plot(c3_est)
title('Estimated values with noise reduction')
legend('c1','c2','c3')

