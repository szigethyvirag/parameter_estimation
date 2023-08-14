clear all, close all, clc;

data = readmatrix('data.xlsx');

y = data(:,1);
yn1 = data(:,2);
yn2 = data(:,3);
u = data(:,4);

N = length(y);

%% estimate theta1 and theta2
Y1 = y(2:end);
X1 = zeros(N-1,2);
X1(:,1) = y(1:end-1);
X1(:,2) = u(1:end-1);

params1 = inv(transpose(X1)*X1)*transpose(X1)*Y1;
theta1 = params1(1);
theta2 = params1(2);

%% estimate theta3 and theta4
Y2 = yn1(2:end);
X2 = zeros(N-1,2);
X2(:,1) = yn1(1:end-1);
X2(:,2) = y(1:end-1);

params2 = inv(transpose(X2)*X2)*transpose(X2)*Y2;
theta3 = params2(1);
theta4 = params2(2);

%% estimate theta5, theta6 and theta7
Y3 = yn2(2:end);
X3 = zeros(N-1,3);
X3(:,1) = yn1(1:end-1);
X3(:,2) = y(1:end-1);
X3(:,3) = yn2(1:end-1);

params3 = inv(transpose(X3)*X3)*transpose(X3)*Y3;
theta5 = params3(1);
theta6 = params3(2);
theta7 = params3(3);

%% simulation:

y_est = zeros(N,1);
yn1_est = zeros(N,1);
yn2_est = zeros(N,1);

for k = 2:N
    y_est(k) = theta1*y_est(k-1) + theta2*u(k-1);
    yn1_est(k) = theta3*yn1_est(k-1) + theta4*y_est(k-1);
    yn2_est(k) = theta5*yn1_est(k-1) + theta6*y_est(k-1) + theta7*yn2_est(k-1);
end

%% plot:
t = 1:1:N;
figure
hold on
plot(t,y,'b'), plot(t,y_est,'r'), legend('Original', 'Estimated'), title('y');
figure
hold on
plot(t,yn1,'b'), plot(t,yn1_est,'r'), legend('Original', 'Estimated'), title('yn1');
figure
hold on
plot(t,yn2,'b'), plot(t,yn2_est,'r'), legend('Original', 'Estimated'), title('yn2');

%% instrumental variable method:
z = zeros(N,1);
zn1 = zeros(N,1);
zn2 = zeros(N,1);

for k = 2:N
    z(k) = theta1*z(k-1) + theta2*u(k-1); 
    zn1(k) = theta3*zn1(k-1) + theta4*z(k-1);
    zn2(k) = theta5*zn1(k-1) + theta6*z(k-1) +theta7*zn2(k-1);
end

xi1 = zeros(N-1,2); 
xi1(:,1) = z(1:end-1);
xi1(:,2) = u(1:end-1);
params1_iv = inv(transpose(xi1)*X1)*transpose(xi1)*Y1; %contains theta1_iv, theta2_iv

xi2 = zeros(N-1,2);
xi2(:,1) = zn1(1:end-1);
xi2(:,2) = z(1:end-1);
params2_iv = inv(transpose(xi2)*X2)*transpose(xi2)*Y2; %contains theta3_iv, theta4_iv

xi3 = zeros(N-1,3);
xi3(:,1) = zn1(1:end-1);
xi3(:,2) = z(1:end-1);
xi3(:,3) = zn2(1:end-1);
params3_iv = inv(transpose(xi3)*X3)*transpose(xi3)*Y3; %contains theta1_i5, theta6_iv, theta7_iv

%% simulation 2:

y_iv_est = zeros(N,1);
yn1_iv_est = zeros(N,1);
yn2_iv_est = zeros(N,1);

for k = 2:N
    y_iv_est(k) = params1_iv(1)*y_iv_est(k-1) + params1_iv(2)*u(k-1);
    yn1_iv_est(k) = params2_iv(1)*yn1_iv_est(k-1) + params2_iv(2)*y_iv_est(k-1);
    yn2_iv_est(k) = params3_iv(1)*yn1_iv_est(k-1) + params3_iv(2)*y_iv_est(k-1) + params3_iv(3)*yn2_iv_est(k-1);
end
%% plot:
t = 1:1:N;
figure
hold on
plot(t,y,'b'), plot(t,y_iv_est,'r'), legend('Original', 'Estimated with iv'), title('y');
figure
hold on
plot(t,yn1,'b'), plot(t,yn1_iv_est,'r'), legend('Original', 'Estimated with iv'), title('yn1');
figure
hold on
plot(t,yn2,'b'), plot(t,yn2_iv_est,'r'), legend('Original', 'Estimated with iv'), title('yn2');

% it doesn't seem to make a difference :(




