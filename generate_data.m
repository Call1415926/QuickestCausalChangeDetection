function output = generate_data(A, mu, Sigma,n, j,value)
%UNTITLED generate causal network data by A mu Sigma
% input: A: weighted matrix A_{ij} neq 0 means Aj is the parent of A_i. p times p
%        mu:  mean of noise.  p times 1
%        Sigma: variance of noise. diagonal p times p
%        n : number of generate data samples
%        j, value: intervention node and intervention value 
%output: n times p matrix. 
p = size(A,1);
I = eye(p);
if j ~= 0
    mu(j) = value;
    Sigma(j,j) = 0;
    A(j,:) = 0;
    
end
mu_X = (I - A) \ mu;
Sigma_X =(I-A) \ Sigma * (inv(I-A))';

output = mvnrnd(mu_X,Sigma_X,n);

end