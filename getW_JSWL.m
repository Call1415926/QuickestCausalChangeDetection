function W = getW_JSWL(tt, p, w, n, tau, b, ...
    A, mu, Sigma, A_oc, mu_oc, Sigma_oc)
%GETW_JSWL Calculate the monitoring statistics of JS-WL-CuSum.
%
% This function implements the fixed-window JS-WL-CuSum using the
% global-mean positive-part James-Stein estimator.
%
% Input:
%   tt        : replication index
%   p         : number of nodes/data streams
%   w         : fixed sliding-window length
%   n         : maximum monitoring length
%   tau       : change time; set tau = Inf for the pre-change condition
%   b         : detection threshold; set b = Inf to return the full path
%   A         : pre-change weighted adjacency matrix
%   mu        : pre-change noise mean
%   Sigma     : pre-change noise covariance matrix
%   A_oc      : post-change weighted adjacency matrix
%   mu_oc     : post-change noise mean
%   Sigma_oc  : post-change noise covariance matrix
%
% Output:
%   W         : 1-by-n JS-WL-CuSum monitoring-statistic sequence
%
% The observations are standardized using the known pre-change
% parameters so that the pre-change distribution is N(0,I).

if size(A, 1) ~= p || size(A, 2) ~= p
    error('The dimension of A must agree with p.');
end

if w >= n
    error('The window length w must be smaller than n.');
end

if w < 1
    error('The window length w must be positive.');
end

if p < 3
    error('The James-Stein estimator requires p >= 3.');
end

W = zeros(1, n);

%% Generate the complete observational sequence

if isinf(tau) || tau > n
    % All observations are generated from the pre-change model
    data = generate_data( ...
        A, mu, Sigma, n, 0, 0);

elseif tau <= 1
    % All observations are generated from the post-change model
    data = generate_data( ...
        A_oc, mu_oc, Sigma_oc, n, 0, 0);

else
    % Observations 1,...,tau-1 are pre-change
    % Observations tau,...,n are post-change
    n_pre = tau - 1;
    n_post = n - n_pre;

    data_pre = generate_data( ...
        A, mu, Sigma, n_pre, 0, 0);

    data_post = generate_data( ...
        A_oc, mu_oc, Sigma_oc, n_post, 0, 0);

    data = [data_pre; data_post];
end

%% Standardize data using the known pre-change parameters

I = eye(p);

% Pre-change observational mean:
%
%   mu_X = (I-A)^{-1} mu.

mu_X = (I - A) \ mu;

% For a column observation X:
%
%   Z = Sigma^{-1/2}(I-A)(X-mu_X).
%
% The following is the equivalent row-observation implementation.

Sigma_inv_sqrt = inv(sqrtm(Sigma));

central_matrix = ...
    (I - A)' * Sigma_inv_sqrt';

Z = (data - mu_X') * central_matrix;

% Under the pre-change model, each row of Z follows N(0,I).

%% Calculate the fixed-window JS-WL-CuSum statistic

for t = (w + 1):n

    % Use only the previous w observations to estimate the unknown
    % post-change mean. The current observation Z(t,:) is not included.
    window_data = Z((t - w):(t - 1), :);

    % Component-wise sample mean vector
    mean_window = mean(window_data, 1);

    % Global mean of all components
    global_mean = mean(mean_window);

    % Projection onto span(1)
    mean_target = global_mean * ones(1, p);

    % Difference from the global-mean shrinkage target
    mean_difference = mean_window - mean_target;

    norm_squared = sum(mean_difference.^2);

    % Global-mean positive-part James-Stein shrinkage factor:
    %
    % a = max{0, 1 - (p-3)/(w*||mean_window-mean_target||^2)}.

    if norm_squared <= eps
        shrinkage_factor = 0;
    else
        shrinkage_factor = max( ...
            0, ...
            1 - (p - 3) / (w * norm_squared));
    end

    % James-Stein estimate of the post-change mean
    theta_hat = ...
        shrinkage_factor * mean_difference ...
        + mean_target;

    %% Gaussian plug-in log-likelihood-ratio increment

    current_data = Z(t, :);

    % For N(theta_hat,I) versus N(0,I):
    %
    % log[f_theta_hat(Z_t)/f_0(Z_t)]
    % = theta_hat' Z_t - 0.5 ||theta_hat||^2.

    statistic = ...
        theta_hat * current_data' ...
        - 0.5 * sum(theta_hat.^2);

    %% WL-CuSum update

    % The original update is:
    %
    % S_t = max(S_{t-1},0) + statistic.

    W(t) = max(W(t - 1), 0) + statistic;

    % Normally inactive when b = Inf
    if W(t) > b
        break;
    end
end

%% Record the replication index

fid = fopen('tt.txt', 'a');
fprintf(fid, '%g\n', tt);
fclose(fid);

end