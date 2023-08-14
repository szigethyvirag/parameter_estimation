function Theta = szivi4_hw2task1(X, Y)
%Dimension of X should be R^(nxp)
%Dimension of Y should be R^n
%Dimension of theta will be R^p
    Theta = inv(transpose(X)*X)*transpose(X)*Y;
end

