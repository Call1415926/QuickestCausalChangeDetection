function [mu_X,Sigma_X] = CalDisX(A, mu, Sigma, j, value)
%calculate distribution of X
%input:  A,mu,Sigma: parameter of DAG
%        j,value: intervention node and intervention value
%output: mu_X, Sigma_X: mean and covariance of X
p = size(A,1);
I = eye(p);
if j ~= 0
    mu(j) = value;
    Sigma(j,j) = 0;
    A(j,:) = 0;
end
mu_X = (I - A) \ mu;
Sigma_X =(I-A) \ Sigma * (inv(I-A))';
end

