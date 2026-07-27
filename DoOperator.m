function [A_do, mu_do,Sigma_do] = DoOperator(A,mu,Sigma,j,value)
%DOOPERATOR: intervention operator to A mu Sigma. j and value are
%intervetino node and intervention value.
if j == 0    
    A_do = A;
    mu_do = mu;
    Sigma_do = Sigma;
else
    A_do = A;
    mu_do = mu;
    Sigma_do = Sigma;
    A_do(j,:) = 0;
    mu_do(j) = value;
    Sigma_do(j,j) = 0;
end


end

