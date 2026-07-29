%% Set parameters for Synthetic Experiment

p = 6;                          % Number of nodes
weightRange = [1, 2];           % Range of nonzero elements of A
maxInDegree = 2;                % Maximum in-degree
muRange = [-1, 1];              % Range of noise means
SigmaRange = [1/2, 2];          % Range of noise variances
ChangeRange = [0.1, 2];         % Range of possible change magnitudes
myChange = 0.1;                 % Actual change magnitude
delta = 1;                      % Intervention gap
w = 80;                         % Window length
q = 40;                         % Number of exploration times
n = 20000;                      % Maximum run length
eta = 1;
Nrep = 1000;                    % Number of replications

[A, mu, Sigma] = generateSparseDAG( ...
    p, weightRange, maxInDegree, muRange, SigmaRange);

%% Randomly select a lower-triangular coefficient to change

[row_set, col_set] = find(tril(true(p), -1));

change_index = randi(length(row_set));

change_row = row_set(change_index);
change_col = col_set(change_index);

% Under a change in A(change_row,change_col), change_col is the
% origin/parent node and is therefore the oracle intervention node.
oraNode = change_col;

A_oc = A;
A_oc(change_row, change_col) = ...
    A(change_row, change_col) + myChange;

mu_oc = mu;
Sigma_oc = Sigma;


%% Set parameters for Ecology Experiment
%
% Comment out the Synthetic Experiment block before activating this block.

% p = 11;
% ChangeRange = [0.1, 2];
% myChange = 0.1;
% delta = 1;
% w = 60;
% q = 30;
% n = 20000;
% eta = 1;
% Nrep = 1000;
%
% load('NECinitial.mat');
%
% maxInDegree = max(sum(A ~= 0, 2));
%
% mu_oc = mu;
% Sigma_oc = Sigma;
% A_oc = A;
%
% % Change 1
% A_oc(11,3) = A(11,3) + myChange;
% oraNode = 3;
%
% % Change 2
% % A_oc(11,9) = A(11,9) + myChange;
% % oraNode = 9;


%% Set parameters for Psychology Experiment
%
% Comment out the Synthetic Experiment block before activating this block.

% p = 5;
% ChangeRange = [0.1, 2];
% myChange = 0.1;
% delta = 1;
% w = 25;
% q = 12;
% n = 20000;
% eta = 1;
% Nrep = 1000;
%
% load('mindInitial.mat');
%
% maxInDegree = max(sum(A ~= 0, 2));
%
% mu_oc = mu;
% Sigma_oc = Sigma;
% A_oc = A;
%
% % Change 1
% A_oc(3,1) = A(3,1) + myChange;
% oraNode = 1;
%
% % Change 2
% % A_oc(5,4) = A(5,4) + myChange;
% % oraNode = 4;


%% Set parameters for Multiple-Change Scenario: Case 1
%
% Comment out the Synthetic Experiment block before activating this block.

% A = [
%     0, 0, 0, 0;
%     1, 0, 0, 0;
%     0, 1, 0, 0;
%     0, 0, 1, 0
% ];
%
% A_oc = [
%     0,   0,    0,   0;
%     1.2, 0,    0,   0;
%     0,   1.14, 0,   0;
%     0,   0,    1.1, 0
% ];
%
% mu = ones(4,1);
% Sigma = eye(4);
%
% mu_oc = mu;
% Sigma_oc = Sigma;
%
% p = 4;
% maxInDegree = max(sum(A ~= 0, 2));
% ChangeRange = [0.1, 2];
% delta = 1;
% w = 40;
% q = 20;
% n = 20000;
% eta = 1;
% Nrep = 1000;
%
% % The joint KL divergence is maximized by intervening on node 1.
% oraNode = 1;


%% Set parameters for Multiple-Change Scenario: Case 2
%
% Comment out the Synthetic Experiment block before activating this block.

% A = [
%     0, 0, 0, 0, 0;
%     1, 0, 0, 0, 0;
%     0, 1, 0, 0, 0;
%     1, 0, 0, 0, 0;
%     0, 0, 0, 1, 0
% ];
%
% A_oc = [
%     0,    0,   0, 0,   0;
%     1.15, 0,   0, 0,   0;
%     0,    1.1, 0, 0,   0;
%     1.15, 0,   0, 0,   0;
%     0,    0,   0, 1.1, 0
% ];
%
% mu = ones(5,1);
% Sigma = eye(5);
%
% mu_oc = mu;
% Sigma_oc = Sigma;
%
% p = 5;
% maxInDegree = max(sum(A ~= 0, 2));
% ChangeRange = [0.1, 2];
% delta = 1;
% w = 60;
% q = 30;
% n = 20000;
% eta = 1;
% Nrep = 1000;
%
% % The joint KL divergence is maximized by intervening on node 1.
% oraNode = 1;


%% Calculate intervention values

doValue = CalDoValue( ...
    A, mu, Sigma, delta, ChangeRange);


%% Experimental settings

ARL_set = [
    100, 200, 400, 800, 1000, ...
    2000, 4000, 6000, 8000, 10000
];

epsilon_set = [
    5, 5, 5, 5, 10, ...
    10, 10, 20, 20, 50
];

method_names = {
    'MAX-AI', ...
    'MAX-RI', ...
    'MAX-NI', ...
    'MULTI-AI', ...
    'MULTI-RI', ...
    'MULTI-NI', ...
    'MULTI-Oracle', ...
    'MAX-Oracle', ...
    'GGM', ...
    'JS-WL-CuSum'
};

num_methods = length(method_names);
num_ARL = length(ARL_set);

EDD_set_10methods = zeros(num_methods, num_ARL);
threshold_set = zeros(num_methods, num_ARL);
achieved_ARL_set = zeros(num_methods, num_ARL);

% All getW-type functions should return the complete statistic path.
b_sim = Inf;

% Prevent an infinite loop when the simulated ARL is discrete.
max_bisection_iter = 1000;


%% Start the parallel pool once

existing_pool = gcp('nocreate');
created_pool = isempty(existing_pool);

if created_pool
    parpool(4);
end


%% Calculate monitoring statistics, thresholds, ARLs and EDDs

for methods = 1:num_methods

    fprintf('Running method %d/%d: %s\n', ...
        methods, num_methods, method_names{methods});

    W_IC = zeros(Nrep, n);
    W_OC = zeros(Nrep, n);

    %% Methods 1-6: original MAX/MULTI methods

    if methods <= 6

        if methods == 1
            monit_type = 'max';
            inter_type = 'AI';

        elseif methods == 2
            monit_type = 'max';
            inter_type = 'RI';

        elseif methods == 3
            monit_type = 'max';
            inter_type = 'NI';

        elseif methods == 4
            monit_type = 'multi';
            inter_type = 'AI';

        elseif methods == 5
            monit_type = 'multi';
            inter_type = 'RI';

        else
            monit_type = 'multi';
            inter_type = 'NI';
        end

        % Pre-change statistic paths
        tau = Inf;

        parfor tt = 1:Nrep
            W_IC(tt,:) = getW( ...
                tt, p, maxInDegree, w, q, n, eta, tau, ...
                inter_type, monit_type, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc, ...
                doValue);
        end

        % Post-change statistic paths
        tau = 1;

        parfor tt = 1:Nrep
            W_OC(tt,:) = getW( ...
                tt, p, maxInDegree, w, q, n, eta, tau, ...
                inter_type, monit_type, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc, ...
                doValue);
        end


    %% Method 7: MULTI-Oracle

    elseif methods == 7

        tau = Inf;

        parfor tt = 1:Nrep
            W_IC(tt,:) = getW_multiOra( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc, ...
                doValue, oraNode);
        end

        tau = 1;

        parfor tt = 1:Nrep
            W_OC(tt,:) = getW_multiOra( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc, ...
                doValue, oraNode);
        end


    %% Method 8: MAX-Oracle

    elseif methods == 8

        tau = Inf;

        parfor tt = 1:Nrep
            W_IC(tt,:) = getW_maxOra( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc, ...
                doValue, oraNode);
        end

        tau = 1;

        parfor tt = 1:Nrep
            W_OC(tt,:) = getW_maxOra( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc, ...
                doValue, oraNode);
        end


    %% Method 9: GGM baseline

    elseif methods == 9

        tau = Inf;

        parfor tt = 1:Nrep
            W_IC(tt,:) = getW_GGM( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc);
        end

        tau = 1;

        parfor tt = 1:Nrep
            W_OC(tt,:) = getW_GGM( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc);
        end


    %% Method 10: JS-WL-CuSum baseline

    elseif methods == 10

        tau = Inf;

        parfor tt = 1:Nrep
            W_IC(tt,:) = getW_JSWL( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc);
        end

        tau = 1;

        parfor tt = 1:Nrep
            W_OC(tt,:) = getW_JSWL( ...
                tt, p, w, n, tau, b_sim, ...
                A, mu, Sigma, ...
                A_oc, mu_oc, Sigma_oc);
        end
    end


    %% Calibrate threshold and calculate EDD

    EDD_set = zeros(1, num_ARL);

    for i = 1:num_ARL

        target_ARL = ARL_set(i);
        epsilon_ARL = epsilon_set(i);

        b_min = 0;
        b_max = max(W_IC(:));

        best_b = b_min;
        best_ARL = 0;
        best_error = Inf;

        for iter = 1:max_bisection_iter

            b_trial = (b_min + b_max) / 2;

            % Calculate the simulated pre-change ARL
            logicalMatrix = W_IC > b_trial;

            [hasAlarm, alarmTime] = ...
                max(logicalMatrix, [], 2);

            alarmTime(~hasAlarm) = size(W_IC,2);

            simulated_ARL = mean(alarmTime);

            current_error = abs( ...
                simulated_ARL - target_ARL);

            % Save the threshold producing the closest ARL
            if current_error < best_error
                best_error = current_error;
                best_b = b_trial;
                best_ARL = simulated_ARL;
            end

            % Stop if the requested tolerance is attained
            if current_error < epsilon_ARL
                break;
            end

            % ARL increases monotonically with the threshold
            if simulated_ARL > target_ARL
                b_max = b_trial;
            else
                b_min = b_trial;
            end

            % Stop if numerical precision prevents further progress
            if abs(b_max - b_min) ...
                    <= 1e-12 * max(1, abs(b_max))
                break;
            end
        end

        b_selected = best_b;

        threshold_set(methods,i) = b_selected;
        achieved_ARL_set(methods,i) = best_ARL;

        % Calculate the post-change EDD under the selected threshold
        logicalMatrix = W_OC > b_selected;

        [hasAlarm, alarmTime] = ...
            max(logicalMatrix, [], 2);

        alarmTime(~hasAlarm) = size(W_OC,2);

        EDD_set(i) = mean(alarmTime);
    end

    EDD_set_10methods(methods,:) = EDD_set;

    fprintf('Completed: %s\n', method_names{methods});
end


%% Close the parallel pool if this script created it

if created_pool
    delete(gcp('nocreate'));
end



%% Plot EDD versus ARL

figure;
hold on;

for i = 1:num_methods
    plot( ...
        ARL_set, ...
        EDD_set_10methods(i,:), ...
        'DisplayName', method_names{i});
end

legend('Location', 'best');
set(gca, 'XScale', 'log');

xlabel('ARL');
ylabel('EDD');

grid on;
hold off;