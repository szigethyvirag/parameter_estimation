load("paramest_hw1_task2.mat");

figure
plot(X(:,1),X(:,2),'.r')
hold on
X_mean = mean(X,1);
plot(X_mean(1),X_mean(2), '*k')

X_cov = cov(X);
z = mvnpdf(X,X_mean,X_cov);

[X_grid,Y_grid] = meshgrid(-2:0.1:6, -2:0.1:6);
Z_grid = griddata(X(:,1),X(:,2),z,X_grid,Y_grid);
mesh(X_grid,Y_grid,Z_grid, 'FaceAlpha', 0.5);
grid on
%rotate3d on
view([36,35])

