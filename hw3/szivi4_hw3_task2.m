clear all
close all

ts = 0.01;

beta = 1.1; % susceptible -> infected 
gamma = 0.4; % recovery rate

N = 1000;
S = zeros(1, N); %susceptibles
I = zeros(1, N); %infected
R = zeros(1, N); %recovered

S(1) = 0.9999;
I(1) = 0.0001;
R(1) = 0;

for k = 2:N
    S(k) = S(k-1) + ts*(-beta*S(k-1)*I(k-1));
    I(k) = I(k-1) + ts*(beta*S(k-1)*I(k-1)-gamma*I(k-1));
    R(k) = R(k-1) + ts*(gamma*I(k-1));
end

figure
hold on
plot(S)
plot(I)
plot(R)
grid on
legend({'susceptible','infected','recovered'},'Location','northwest')
%% 

load("sir.csv")
ts = 0.01;

S_measured = sir(:,1);
I_measured = sir(:,2);
R_measured = sir(:,3);

X = zeros(500,2);

X(:,1) = S_measured(1:end-1);
X(:,2) = -ts * (S_measured(1:end-1) .* I_measured(1:end-1));
Y = S_measured(2:end);
beta_pred = (X' * X) \ (X' * Y);

X(:,1) = R_measured(1:end-1);
X(:,2) = ts * I_measured(1:end-1);
Y = R_measured(2:end);
gamma_pred = (X' * X) \ (X' * Y);

disp(beta_pred);
disp(gamma_pred);
