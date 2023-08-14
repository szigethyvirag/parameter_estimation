clear all
close all

load("paramest_hw1_task1.mat");

y = zeros(size(u));
e = normrnd(0,1.2,size(u));
N = 100;
x = 1:100;
p = [1.5; -0.7; -0.25; 1.64];
ep = [-1.6; 0.75];

%first element:
y(1) = p(1)*0+p(2)*0+p(3)*0+p(4)*0+ep(1)*0+ep(2)*0;
%second element:
y(2) = p(1)*y(1)+p(2)*0+p(3)*u(1)+p(4)*0+ep(1)*e(1)+ep(2)*0;
%all the other elements:
for k = 3:N
    y(k) = p(1)*y(k-1)+p(2)*y(k-2)+p(3)*u(k-1)+p(4)*u(k-2)+ep(1)*e(k-1)+ep(2)*e(k-2);
end


%parameter estimation with LSQ:
X = [[0; 0; y(2:99)],[0; 0; y(1:98)],[0; 0; u(2:99)],[0; 0; u(1:98)]];
theta = (X'*X)\(X'*y);
disp('distances: ' )
disp(p-theta)

%signal with the estimated parameters:
y_est = zeros(size(u));

%first element:
y_est(1) = theta(1)*0+theta(2)*0+theta(3)*0+theta(4)*0;
%second element:
y_est(2) = theta(1)*y_est(1)+theta(2)*0+theta(3)*u(1)+theta(4)*0;
%all the other elements:
for k = 3:N
    y_est(k) = theta(1)*y_est(k-1)+theta(2)*y_est(k-2)+theta(3)*u(k-1)+theta(4)*u(k-2);
end

figure
hold on
plot(x,y)
plot(x,y_est)
grid on
legend('original ARMAX','Signal with the estimated parameters')


