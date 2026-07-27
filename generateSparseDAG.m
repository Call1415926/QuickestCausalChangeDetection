function [A,mu,Sigma] = generateSparseDAG(p, weightRange, maxInDegree, muRange, SigmaRange)
%%%generate sparse DAG A with lower triangle matrix, mu ,diagonal Sigma.
%%%input: p: number of nodes
%         weightRange: weight range of A
%         muRange: range of mu [min,max]
%         SigmaRange: range of sigma [min,max]
%         maxInDegree: d
%output: A: p times p lower trianle matrix
%        mu: p times 1 mean vector
%        Sigma: p times p diagnal matrix

    A = zeros(p);  
    mu = zeros(p,1);
    Sigma = zeros(p);

    for i = 2:p
        numEdges = randi([1, min(maxInDegree, i-1)]);  
        predecessors = randperm(i-1, numEdges);  

        for j = 1:numEdges
            A(i, predecessors(j)) = rand() * (weightRange(2) - weightRange(1)) + weightRange(1);
        end
    end
    for i = 1:p
        mu(i) = rand() * (muRange(2) - muRange(1)) + muRange(1);
        Sigma(i,i) = rand() * (SigmaRange(2) - SigmaRange(1)) + SigmaRange(1);
    end


end


