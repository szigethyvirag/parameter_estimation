close all, clear all;

N = 1000;
sigma1 = 0.6;
sigma2 = 0.8;
e = [normrnd(0,sigma1,[1,N/2-1]) , normrnd(0,sigma2,[1,N/2+1])];

theta1 = 0.3;
theta2 = 1.2;

u = normrnd(0,1.4,[1,N]); % input
y = zeros([1,N]); % output

for k = 2:N
    y(k) = theta1*y(k-1) + theta2*u(k-1)+e(k);
end

log_likelihood(u,y,theta1,theta2,sigma1,sigma2)

%%
theta1_values = 0.2:0.2:1.8;
theta2_values = 0.2:0.2:1.8;
[xx, yy] = meshgrid(theta1_values, theta2_values);

zz = [];
for j = theta2_values
    z_tmp = [];
    for k = theta1_values
        z_tmp = [z_tmp, log_likelihood(u,y,k,j,sigma1,sigma2)];
    end
    zz = [zz;z_tmp];
end

figure 
grid on
surf(xx,yy,zz)
title('log-likelihood surface')
xlabel('theta1')
ylabel('theta2')
zlabel('log-likelihood')
%%
phat = mle(y)
%%
N = 1000;
sigma1 = 0.6;
sigma2 = 0.8;

theta1 = 0.3;
theta2 = 1.2;

u = normrnd(0,1.4,[1,N]); % input
y_c = zeros([1,N]); % 1x1000

theta_estimated_values_maxl = [];
theta_estimated_values_leastsq = [];

for j = 1:100
    e = [normrnd(0,sigma1,[1,N/2-1]) , normrnd(0,sigma2,[1,N/2+1])];
    for k = 2:N
        y_c(k) = theta1*y_c(k-1) + theta2*u(k-1)+e(k);
    end
    x_tmp = [y_c;u];
    theta_estimated_values_maxl = [theta_estimated_values_maxl; mle(y_c)];
    est_theta = szivi4_hw2_task1(transpose(x_tmp), transpose(y)); %3x1000 1x1000
    theta_estimated_values_leastsq = [theta_estimated_values_leastsq; transpose(est_theta)]; %1x1000 
end

cov_maxl = cov(theta_estimated_values_maxl)
cov_leastsq = cov(theta_estimated_values_leastsq)

%%

function fy = log_likelihood(u, y, theta1, theta2, sigma1, sigma2) 
    fy = log(1/(sigma1*sqrt(2*pi))) + log(exp(-(y(1)-theta1*0-theta2*0)^2/(2*sigma1^2)));
    for k = 2:499
        fy = fy + log(1/(sigma1*sqrt(2*pi))) + log(exp(-(y(k)-theta1*y(k-1)-theta2*u(k-1))^2/(2*sigma1^2)));
    end
    for k = 500:1000
        fy = fy + log(1/(sigma2*sqrt(2*pi))) + log(exp(-(y(k)-theta1*y(k-1)-theta2*u(k-1))^2/(2*sigma2^2)));
    end
end