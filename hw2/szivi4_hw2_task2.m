%loss function depends on the amplitude of noise:
clear all

loss_function = [];
sigma_values = 0.2:0.2:0.8;

for sigma = sigma_values
    
    %data generation:
    n = 15;
    x = 10.*rand(n,2)+1;
    x_tmp = [ones(n,1),x];
    theta = 10.*rand(1,3);
    %sigma = 1;
    e = normrnd(0, sigma, [1,n]);
    y = x_tmp*transpose(theta) + transpose(e);

    %3d plot:
    figure
    plot3(x(:,1),x(:,2), y,'o')
    title(['Data with sigma ', num2str(sigma)])
    xlabel('x1')
    ylabel('x2')
    zlabel('y')
    grid on

    %Least square estimation:
    est_theta = szivi4_hw2_task1(x_tmp, y); %40x3 40x1 -> 3

    %Loss function:
    Loss = transpose(y - x_tmp*est_theta)*(y - x_tmp*est_theta);
    loss_function = [loss_function, Loss];
end

figure 
plot(sigma_values, loss_function,'*-')
grid on
title('Loss function with different sigma values')
xlabel('sigma')
ylabel('loss')

%%
%loss function depends on the sample size:
clear all

loss_function = [];
n_values = 10:10:40;

for n = n_values
    
    %data generation:
    %n = 15;
    x = 10.*rand(n,2)+1;
    x_tmp = [ones(n,1),x];
    theta = 10.*rand(1,3);
    sigma = 0.2;
    e = normrnd(0, sigma, [1,n]);
    y = x_tmp*transpose(theta) + transpose(e);

    %3d plot:
    figure
    plot3(x(:,1),x(:,2), y,'o')
    title(['Data with n ', num2str(n)])
    xlabel('x1')
    ylabel('x2')
    zlabel('y')
    grid on

    %Least square estimation:
    est_theta = szivi4_hw2_task1(x_tmp, y);

    %Loss function:
    Loss = transpose(y - x_tmp*est_theta)*(y - x_tmp*est_theta);
    loss_function = [loss_function, Loss];
end

figure 
plot(n_values, loss_function,'*-')
grid on
title('Loss function with different n values')
xlabel('n')
ylabel('loss')

%%
%all the upper combined togeter:
%loss function depends on the amplitude of noise and the the sample size:
clear all

n_values = 10:10:40;
sigma_values = 0.2:0.2:0.8;
zz = [];

for n = n_values
    loss_function = [];
    for sigma = sigma_values
    
        %data generation:
        %n = 15;
        x = 10.*rand(n,2)+1;
        x_tmp = [ones(n,1),x];
        theta = 10.*rand(1,3);
        %sigma = 0.2;
        e = normrnd(0, sigma, [1,n]);
        y = x_tmp*transpose(theta) + transpose(e);

        %Least square estimation:
        est_theta = szivi4_hw2_task1(x_tmp, y);

        %Loss function:
        Loss = transpose(y - x_tmp*est_theta)*(y - x_tmp*est_theta);
        loss_function = [loss_function, Loss];
    end
    zz = [zz; loss_function];
end

[xx, yy] = meshgrid(sigma_values, n_values);
figure 
grid on
surf(xx,yy,zz)
title({['dependency of the loss function'] ['on the amplitude of noise and the the sample size']})
xlabel('sigma')
ylabel('n')
zlabel('loss')
