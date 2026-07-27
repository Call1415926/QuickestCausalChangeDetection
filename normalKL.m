function  distance = normalKL(mu,Sigma)
%NORMALKL calculate the KL divergence with the stand norm distribution.
p = length(mu);

distance = sum(mu.^2) + trace(Sigma) - p - log(det(Sigma));
distance = distance/2;

end

