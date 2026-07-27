function doValue = CalDoValue(A,mu,Sigma,delta,ChangeRange)
%CALDOVALUE calculate intervetnion value
%input: A,mu,Sigma: parameter of DAG
%       delta: intervention gap
%       ChaneRange: [\Delta_{min}, \Delta_{max}]
%output: intervetion value.

p = size(A,1);
doValue = zeros(p,1);
I = eye(p);
E = inv(I - A);

all_set = 1:p;
for j = 1:p   
anc_set = find(E(j,1:(j-1))); 
no_anc_set = setdiff(all_set, [anc_set,j]);
if isempty(no_anc_set)  
    continue
end
    
        for l = [0, anc_set] 
            if l == 0   
                [A_do, mu_do,Sigma_do] = DoOperator(A,mu,Sigma, l);
            else
                [A_do, mu_do,Sigma_do] = DoOperator(A,mu,Sigma, l, doValue(l));
            end

            E_do = inv(I - A_do);
            part1 = (E_do(j,:) * mu_do)^2;
            part2 = E_do(j,:).^2 * diag(Sigma_do);
            temp = diag(Sigma);
            part3 = 2*delta * max(temp(no_anc_set))/ChangeRange(1)^2;
            doValue(j) = max(doValue(j), sqrt(part1+part2+part3));

        end
   

end


end

