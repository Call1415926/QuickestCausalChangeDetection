%% Set parameter for Synthetic Experiment
p = 6;  %  number of nodes
weightRange = [1, 2];  %   range of A
maxInDegree = 2;  % d
muRange = [-1,1]; %%range of mu
SigmaRange = [1/2, 2]; %%range of Sigma
ChangeRange = [0.1, 2]; %% range of change
myChange = 0.1; %% Delta
delta = 1; %%% delta
w = 80; 
q = 40;  
n = 20000;   %% max runlength 
eta = 1; 
b = Inf; 
Nrep = 1000;  %%%%numner of replications
[A,mu,Sigma] = generateSparseDAG(p, weightRange, maxInDegree,muRange,SigmaRange); %%%generate A mu Sigma

%%%%%%random choose a edge to change 
[row,col] = find(tril(true(p),-1));

idx = randi(length(row));
row = row(idx);
col = col(idx);

A_oc = A;
A_oc(row,col) = A(row,col) + myChange;  %%%set post parameter
mu_oc = mu;
Sigma_oc = Sigma;

%% Set parameter for ecology experiment
% p = 11;  %  number of nodes
% ChangeRange = [0.1, 2]; %% range of change
% myChange = 0.1; %% Delta
% delta = 1; %%% delta
% w = 60; %%%% 
% q = 30; %%% 
% n = 20000;   %% max runlength 
% eta = 1; %%%
% b = Inf; 
% Nrep = 1000;  %%%%numner of replications
% load('NECinitial.mat')  %%% pre-change parameter setting
% 
% mu_oc = mu;
% Sigma_oc = Sigma;
% A_oc = A;
% A_oc(11,3) = A(11,3) + myChange;  %%%change 1
% %A_oc(11,9) = A(11,9) + myChange;  %%%change 2

%% Set parameter for psychology experiment
% p = 5;  %  number of nodes
% ChangeRange = [0.1, 2]; %% range of change
% myChange = 0.1; %% Delta
% delta = 1; %%% delta
% w = 25; 
% q = 12; 
% n = 20000;   %% max runlength 
% eta = 1; 
% b = Inf; 
% Nrep = 1000;  %%%%numner of replications
% load('mindInitial.mat')  %%% pre-change parameter setting
% 
% mu_oc = mu;
% Sigma_oc = Sigma;
% A_oc = A;
% A_oc(3,1) = A(3,1) + myChange;  %%%change 1
% %A_oc(5,4) = A(5,4) + myChange;  %%%change 2


%% Set parmater for multiple changes secnario
%%% Case 1
%A = [0,0,0,0;
%    1,0,0,0;
%    0,1,0,0;
%    0,0,1,0];
%A_oc = [0,0,0,0;
%    1.2,0,0,0;
%    0,1.14,0,0;
%    0,0,1.1,0];
%mu = ones(4,1);
%Sigma = eye(4);
% p = 4;  %  number of nodes
% ChangeRange = [0.1, 2]; %% range of change
% delta = 1; %%% delta
% w = 40; 
% q = 20;  
% n = 20000;   %% max runlength 
% eta = 1; 
% b = Inf; 
% Nrep = 1000;  %%%%numner of replications

%%% Case 2
% A = [0,0,0,0,0;
%     1,0,0,0,0;
%     0,1,0,0,0;
%     1,0,0,0,0;
%     0,0,0,1,0];
% A_oc = [0,0,0,0,0;
%     1.15,0,0,0,0;
%     0,1.1,0,0,0;
%     1.15,0,0,0,0;
%     0,0,0,1.1,0];
% mu = ones(5,1);
% Sigma = eye(5);
% mu_oc = mu;
% Sigma_oc = Sigma;
% p = 5;  %  number of nodes
% ChangeRange = [0.1, 2]; %% range of change
% delta = 1; %%% delta
% w = 60;  
% q = 30;  
% n = 20000;   %% max runlength 
% eta = 1; 
% b = Inf; 
% Nrep = 1000;  %%%%numner of replications


%% Calculate intervention value
doValue = CalDoValue(A,mu,Sigma,delta,ChangeRange);

%% Calculate monitoring statistics, ARL, EDD for four methods.

ARL_set = [100,200,400,800,1000,2000,4000,6000,8000,10000]; %%% ARL set. w  must be smaller than min_ARL, or the ARL need to be reset
epsilon_set = [5,5,5,5,10,10,10,20,20,50];
EDD_set_4methods = zeros(4,length(ARL_set));
for methods = 1:6
    if methods == 1
        %%% MAX-AI
        monit_type = 'max';
        inter_type = 'AI';
    end
    
    if methods == 2
        %%% MAX-RI
        monit_type = 'max';
        inter_type = 'RI';
    end
    
    if methods == 3
        %%% MAX-NI
        monit_type = 'max';
        inter_type = 'NI';
    end
    
    if methods == 4
        %%% Multi-AI
        monit_type = 'multi';
        inter_type = 'AI';
    end

    if methods == 5
        %%% multi-RI
        monit_type = 'multi';
        inter_type = 'RI';
    end

    if methods == 6
        %%% multi-NI
        monit_type = 'multi';
        inter_type = 'NI';
    end
    
    tau = Inf;  %%%% calculate monitoring statistics for pre-change condition
    W_IC = zeros(Nrep,n);
    parpool(4);
    parfor tt = 1:Nrep
       W_IC(tt,:)  = getW(tt, p,maxInDegree, w,q,n,eta,tau, inter_type, monit_type,b,A,mu,Sigma,A_oc,mu_oc,Sigma_oc,doValue);
    end
    delete(gcp('nocreat'));
    
    tau = 1; %%%% calculate monitoring statistics for post-change condition
    W_OC = zeros(Nrep,n);
    parpool(4);
    parfor tt = 1:Nrep
       W_OC(tt,:)  = getW(tt, p,maxInDegree, w,q,n,eta,tau, inter_type, monit_type,b,A,mu,Sigma,A_oc,mu_oc,Sigma_oc,doValue);
    end
    delete(gcp('nocreat'));
    
    %%%calculate threshold b and corresponding EDD
    EDD_set = zeros(1,length(ARL_set));
    for i = 1:length(ARL_set)
        ARL = ARL_set(i);
        epsilon = epsilon_set(i);
        %%%%pre-change
        b_max = max(W_IC(:));
        b_min = 0;
        while true
            b = (b_min + b_max)/2;
            logicalMatrix = W_IC > b;
            [idx, col] = max(logicalMatrix, [], 2);
            col(~idx) = size(W_IC,2);  
            myARL = mean(col);
            if abs(myARL - ARL) < epsilon
                break
            end
            if myARL > ARL
                b_max = b;
            else
                b_min = b;
            end
        end
        %%%%post-change 
        logicalMatrix = W_OC > b;
        [idx, col] = max(logicalMatrix, [], 2);
        col(~idx) = size(W_OC,2);  
        EDD = mean(col);
        EDD_set(i) = EDD;
    end
   EDD_set_4methods(methods,:) = EDD_set;
end

%% plot EDD v.s. ARL

figure;
hold on;  
names = {'MAX-AI','MAX-RI','MAX-NI','MULTI-AI','MULTI-RI','MULTI-NI'}; %% legend name
for i = 1:6
    plot(ARL_set, EDD_set_4methods(i,:), 'DisplayName', names{i}); 
end
legend;
set(gca, 'XScale', 'log'); %% log transform
xlabel('ARL');
ylabel('EDD');
hold off;

