clear all, close all, clc;

in = readmatrix('pharmacokinetics.xlsx');
x1 = in(:,1);
x2 = in(:,2);
x3 = in(:,3);
x4 = in(:,4);

vm = 10;
x1_0 = 0.27;
x2_0 = 0.2;
x3_0 = 0.2;
x4_0 = 0.2;
ts = 0.1;

%% estimate alpha2 and beta2 

N = length(x1);
Y2 = zeros(N,1);
Y4 = zeros(N,1);

Y2(1) = (x2(1) - x2_0)/ts;
Y4(1) = (x4(1) - x4_0)/ts;
Y2(2:end) = (x2(2:end) - x2(1:end-1))/ts;
Y4(2:end) = (x4(2:end) - x4(1:end-1))/ts;

X2 = x1-x2;
X4 = x3-x4;

alpha2 = inv(transpose(X2)*X2)*transpose(X2)*Y2
beta2 = inv(transpose(X4)*X4)*transpose(X4)*Y4

%% estimate alpha1 and beta1

x1_tmp = zeros(N,1); % for x1(k-1)
x3_tmp = zeros(N,1); % for x3(k-1)
x1_tmp(1) = x1_0;
x3_tmp(1) = x3_0;
x1_tmp(2:end) = x1(1:end-1);
x3_tmp(2:end) = x3(1:end-1);

% from derivation on paper constructing X:
X = zeros(2*N,10);
X(1:N,1) = x1 - x1_tmp;                 % kc 
X(1:N,2) = 0.1*x1.*x1 - 0.1*x1.*x2;     % alpha1
X(1:N,3) = x1.*x3 - x1_tmp.*x3;         % kc/ka
X(1:N,4) = 0.1*x1 - 0.1*x2;             % alpha1*kc
X(1:N,5) = 0.1*x1 - 0.1*x2.*x3;         % alpha1*kc/ka
X(N+1:end,6) = x3 - x3_tmp;             % ka
X(N+1:end,7) = 0.1*x3.*x3 - 0.1*x3.*x4; % beta1
X(N+1:end,8) = x1.*x3 - x1.*x3_tmp;     % ka/kc
X(N+1:end,9) = 0.1*x3 - 0.1*x4;         % beta1*ka
X(N+1:end,10) = 0.1*x1.*x3 - 0.1*x1.*x4;% beta1*ka/kc

Y = zeros(2*N,1);
Y(1:N) = x1.*x1_tmp - x1.*x1 - x1;      % from eq1
Y(N+1:end) = x3.*x3_tmp - x3.*x3 - x3;  % from eq2

% setting the constraint:
D = 4.5;
C = [1 0 0 0 0 1 0 0 0 0];

% two helping matrices
A = zeros(11,11);
A(1:10,1:10) = 2*transpose(X)*X;
A(1:10,11) = transpose(C);
A(11,1:10) = C;

B = zeros(11,1);
B(1:10) = 2*transpose(X)*Y;
B(11) = D;

tetha = A * B; % actually it also includes lambda
kc = tetha(1)
alpha1 = tetha(2)
ka = tetha(6)
beta1 = tetha(7)

%% squared loss

% simulating the system:
x1_est = zeros(N,1);
x2_est = zeros(N,1);
x3_est = zeros(N,1);
x4_est = zeros(N,1); 

% from the first 4 differential equations
x1_est(1) = (alpha1*(x2_0-x1_0) - ka*vm*x1_0/(kc*ka+kc*x3_0+ka*x1_0))*ts + x1_0;
x2_est(1) = (alpha2*(x1_0-x2_0))*ts + x2_0;
x3_est(1) = (beta1*(x4_0-x3_0) - kc*vm*x3_0/(kc*ka+kc*x3_0+ka*x1_0))*ts + x3_0;
x4_est(1) = (beta2*(x3_0-x4_0))*ts + x4_0;

for k = 2:N
    x1_est(k) = (alpha1*(x2_est(k-1)-x1_est(k-1)) - ka*vm*x1_est(k-1)/(kc*ka+kc*x3_est(k-1)+ka*x1_est(k-1)))*ts + x1_est(k-1);
    x2_est(k) = (alpha2*(x1_est(k-1)-x2_est(k-1)))*ts + x2_est(k-1);
    x3_est(k) = (beta1*(x4_est(k-1)-x3_est(k-1)) - kc*vm*x3_est(k-1)/(kc*ka+kc*x3_est(k-1)+ka*x1_est(k-1)))*ts + x3_est(k-1);
    x4_est(k) = (beta2*(x3_est(k-1)-x4_est(k-1)))*ts + x4_est(k-1);
end

x1_sl = 1/N*sum((x1-x1_est).^2)
x2_sl = 1/N*sum((x2-x2_est).^2)
x3_sl = 1/N*sum((x3-x3_est).^2)
x4_sl = 1/N*sum((x4-x4_est).^2)

%%
%plot

t = 0:ts:(N-1)*ts;
figure
hold on
plot(t,x1,'b'), plot(t,x1_est,'r'), legend('Original', 'Estimated'), title('x1');
figure
hold on
plot(t,x2,'b'), plot(t,x2_est,'r'), legend('Original', 'Estimated'), title('x2');
figure
hold on
plot(t,x3,'b'), plot(t,x3_est,'r'), legend('Original', 'Estimated'), title('x3');
figure
hold on
plot(t,x4,'b'), plot(t,x4_est,'r'), legend('Original', 'Estimated'), title('x4');

