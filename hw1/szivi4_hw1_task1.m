load("paramest_hw1_task1.mat");

y = zeros(size(u));
%e = randn(size(u)); %should I do smg else with mean and deviation?
e = normrnd(0,1.2,size(u));
N = 100;
x = 1:100;

%first element:
y(1) = 1.5*0-0.7*0-0.25*0+1.64*0-1.6*0+0.75*0;
%second element:
y(2) = 1.5*y(1)-0.7*0-0.25*u(1)+1.64*0-1.6*e(1)+0.75*0;
%all the other elements:
for k = 3:N
    y(k) = 1.5*y(k-1)-0.7*y(k-2)-0.25*u(k-1)+1.64*u(k-2)-1.6*e(k-1)+0.75*e(k-2);
end

figure
plot(x,y)
grid on

K = 5;
y_transf = movmean(y,K);

hold on 
plot(x,y_transf)

xx = 1:.05:100;
s = spline(x,y,xx);
plot(xx,s)
legend('original ARMAX','moving average','cubic spline')

%derivatives:
d1 = diff(s,1);
d2 = diff(s,2);
d3 = diff(s,3);

figure
grid on
hold on
plot(xx,s)
plot(1:0.05:99.95,d1)
plot(1:0.05:99.9,d2)
plot(1:0.05:99.85,d3)
legend('cubic spline','1st derivative','2nd derivative','3rd derivative')
