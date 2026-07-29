function W = getW_maxOra(tt, p, w, n, tau, b, A, mu, Sigma, A_oc, mu_oc, Sigma_oc, doValue, oraNode)
%GETW_MAXORA Calculate the monitoring statistics of MAX-Oracle.
%
% Input:
%   tt        : replication index
%   p         : number of nodes
%   w         : sliding-window length
%   n         : maximum monitoring length
%   tau       : change time; set tau = Inf for the pre-change condition
%   b         : detection threshold; set b = Inf to return the full path
%   A         : pre-change weighted adjacency matrix
%   mu        : pre-change noise mean
%   Sigma     : pre-change noise covariance matrix
%   A_oc      : post-change weighted adjacency matrix
%   mu_oc     : post-change noise mean
%   Sigma_oc  : post-change noise covariance matrix
%   doValue   : p-dimensional vector of intervention values
%   oraNode   : oracle intervention node
%
% Output:
%   W         : 1-by-n MAX-Oracle monitoring-statistic sequence

if oraNode < 0 || oraNode > p || floor(oraNode) ~= oraNode
    error('oraNode must be an integer in {0,1,...,p}.');
end

if w >= n
    error('The window length w must be smaller than n.');
end

%% Set the fixed oracle intervention

if oraNode == 0
    value = 0;
else
    value = doValue(oraNode);
end

%% Generate all observations under the same oracle intervention

if isinf(tau) || tau > n
    % All observations are pre-change
    data = generate_data( ...
        A, mu, Sigma, n, oraNode, value);

elseif tau <= 1
    % All observations are post-change
    data = generate_data( ...
        A_oc, mu_oc, Sigma_oc, n, oraNode, value);

else
    % Observations before tau are pre-change
    n_pre = tau - 1;
    n_post = n - n_pre;

    data_pre = generate_data( ...
        A, mu, Sigma, n_pre, oraNode, value);

    data_post = generate_data( ...
        A_oc, mu_oc, Sigma_oc, n_post, oraNode, value);

    data = [data_pre; data_post];
end

%% Centralize all observations using the pre-change parameters

if oraNode == 0
    [mu_X, ~] = CalDisX(A, mu, Sigma, 0, 0);

    central_matrix = (eye(p) - A)' ...
        * inv(sqrtm(Sigma))';

    Y = (data - mu_X') * central_matrix;

else
    [mu_X, ~] = CalDisX( ...
        A, mu, Sigma, oraNode, value);

    [A_do, ~, ~] = DoOperator( ...
        A, mu, Sigma, oraNode, value);

    central_matrix = (eye(p) - A_do)' ...
        * inv(sqrtm(Sigma))';

    Y = (data - mu_X') * central_matrix;
end

%% Initialize the node-wise MAX statistics

W_node = zeros(p, n);

%% Calculate the MAX-Oracle monitoring statistic

for t = (w + 1):n

    % The distribution parameters are estimated from the previous
    % w observations, all obtained under the oracle intervention
    window_data = Y((t - w):(t - 1), :);

    mean_hat = mean(window_data, 1);
    cov_hat = cov(window_data, 1);

    current_data = Y(t, :);

    static = zeros(p, 1);

    for l = 1:p

        % The intervened node is deterministic and is not monitored
        if oraNode ~= 0 && l == oraNode
            continue;
        end

        variance_hat = cov_hat(l, l);

        if variance_hat > 0
            static(l) = ...
                log(normpdf( ...
                    current_data(l), ...
                    mean_hat(l), ...
                    variance_hat)) ...
                - log(normpdf(current_data(l), 0, 1));
        end
    end

    % Update the p node-wise CUSUM statistics
    for l = 1:p
        W_node(l, t) = ...
            max(W_node(l, t - 1), 0) + static(l);
    end

    if max(W_node(:, t)) > b
        break;
    end
end

%% Return the maximum node-wise statistic

W = max(W_node, [], 1);

%% Record the replication index

fid = fopen('tt.txt', 'a');
fprintf(fid, '%g\n', tt);
fclose(fid);

end