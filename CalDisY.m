function [mu_Y,Sigma_Y] = CalDisY(A, mu, Sigma, A_oc,mu_oc, Sigma_oc, j, value)
%calculate distribution of Y
%input:  A,mu,Sigma: pre-cahnge parameter of DAG, A_oc mu_oc,Sigma_oc post-change
%parameter of DAG
%        j,value: intervention node and intervention value
%output: mu_X, Sigma_X: mean and covariance of X

p = size(A,1);
I = eye(p);

[mu_X_oc, Sigma_X_oc] = CalDisX(A_oc, mu_oc, Sigma_oc, j, value);
[mu_X, Sigma_X] = CalDisX(A, mu, Sigma, j, value);

if j == 0
   mu_Y = inv(sqrtm(Sigma))*(I-A) * (mu_X_oc - mu_X);
   Sigma_Y = inv(sqrtm(Sigma))*(I-A) * Sigma_X_oc * (I-A)' * inv(sqrtm(Sigma))';
else

    [A,mu,Sigma] = DoOperator(A,mu,Sigma,j,value);
    temp1 = (I-A) * (mu_X_oc - mu_X);
    temp1(j) = [];
    temp2 = Sigma;
    temp2(j,:) = [];
    temp2(:,j) = [];
    temp2 = inv(sqrtm(temp2));
    
    mu_Y = temp2 * temp1;
    
    temp3 = (I-A) * Sigma_X_oc * (I-A)';
    temp3(j,:) = [];
    temp3(:,j) = [];
    
    Sigma_Y = temp2 *temp3 *temp2';
end




end

