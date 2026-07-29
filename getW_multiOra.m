function W = getW_multiOra(tt, p, w, n, tau, b, ...
    A, mu, Sigma, A_oc, mu_oc, Sigma_oc, doValue, oraNode)
%GETW_MULTIORA Calculate the monitoring statistics of MULTI-Oracle.
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
%   W         : 1-by-n MULTI-Oracle monitoring-statistic sequence

if oraNode < 0 || oraNode > p || floor(oraNode) ~= oraNode
    error('oraNode must be an integer in {0,1,...,p}.');
end

if w >= n
    error('The window length w must be smaller than n.');
end

W = zeros(1, n);

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

    % Remove the intervened node
    Y(:, oraNode) = [];
end

multi_dim = size(Y, 2);

%% Calculate the MULTI-Oracle monitoring statistic

for t = (w + 1):n

    % Estimate the post-change mean and covariance using the previous
    % w observations, all obtained under the oracle intervention
    window_data = Y((t - w):(t - 1), :);

    mean_hat = mean(window_data, 1);
    cov_hat = cov(window_data, 1);

    % Check whether the estimated covariance matrix is positive definite
    [~, nonpositive] = chol(cov_hat);

    if nonpositive == 0
        current_data = Y(t, :);

        statistic = ...
            log(mvnpdf(current_data, mean_hat, cov_hat)) ...
            - log(mvnpdf( ...
                current_data, ...
                zeros(1, multi_dim), ...
                eye(multi_dim)));
    else
        statistic = 0;
    end

    % MULTI-type CUSUM update
    W(t) = max(W(t - 1), 0) + statistic;

    if W(t) > b
        break;
    end
end

%% Record the replication index

fid = fopen('tt.txt', 'a');
fprintf(fid, '%g\n', tt);
fclose(fid);

end
