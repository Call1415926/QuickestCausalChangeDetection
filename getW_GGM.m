function W = getW_GGM(tt, p, w, n, tau, b, ...
    A, mu, Sigma, A_oc, mu_oc, Sigma_oc)
%GETW_GGM Calculate the monitoring statistics of the GGM baseline.
%
% This function implements the known pre-change precision-matrix
% version of the sequential GGM change detector proposed by
% Keshavarz et al. (JMLR, 2020).
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
%
% Output:
%   W         : 1-by-n GGM monitoring-statistic sequence
%
% The method uses observational data only. The known pre-change
% observational mean is removed before applying the GGM detector.

if size(A, 1) ~= p || size(A, 2) ~= p
    error('The dimension of A must agree with p.');
end

if w >= n
    error('The window length w must be smaller than n.');
end

if w < 1
    error('The window length w must be positive.');
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

%% Calculate the known pre-change observational mean

% Under the pre-change SEM:
%
%   X = A X + U,   U ~ N(mu, Sigma),
%
% the observational mean is
%
%   mu_X = (I-A)^{-1} mu.

I = eye(p);
mu_X = (I - A) \ mu;

% Center all observations using the known pre-change mean
data_centered = data - mu_X';

%% Calculate the known pre-change precision matrix

% The pre-change observational covariance matrix is
%
%   Sigma_X = (I-A)^{-1} Sigma (I-A)^{-T}.
%
% Therefore, its precision matrix is
%
%   Omega_0 = (I-A)' Sigma^{-1} (I-A).

Omega_0 = (I - A)' * (Sigma \ (I - A));

% Remove small numerical asymmetry
Omega_0 = (Omega_0 + Omega_0') / 2;

if any(diag(Omega_0) <= 0)
    error('The diagonal elements of Omega_0 must be positive.');
end

%% Calculate the standardized precision matrix R

omega_diag = diag(Omega_0);

R_denominator = sqrt(omega_diag * omega_diag');

R = Omega_0 ./ R_denominator;

% The paper uses h_w(r) approximately equal to r^4.
% Thus:
%
% sqrt(sum_{s1,s2} h_w(R_{s1,s2}))
% approximately equals
% sqrt(sum_{s1,s2} R_{s1,s2}^4).

R_norm = sqrt(sum(R(:).^4));

%% Calculate g1(w) and g2(w)

% psi(w/2) is the digamma function.
g1 = log(w / 2) - psi(w / 2);

% psi(1,w/2) is the trigamma function.
g2_squared = psi(1, w / 2) - 2 / w;

if g2_squared <= 0
    error('The calculated value of g2(w)^2 is nonpositive.');
end

g2 = sqrt(g2_squared);

standardization = g2 * R_norm;

if standardization <= 0
    error('The standardization term must be positive.');
end

%% Calculate the sliding-window GGM statistics

% At monitoring time u, use observations
%
%   X_{u-w+1},...,X_u
%
% to calculate the paper's statistic T_{u-w}. The first w entries
% of W remain zero.

for u = (w + 1):n

    window_data = data_centered((u - w + 1):u, :);

    % For every observation x, calculate x' * Omega_0.
    % The s-th column is <x, Omega_{0,:,s}>.
    conditional_score = window_data * Omega_0;

    % Calculate Y_s^{(t,w)} for all s:
    %
    % Y_s = sum_r <X_r,Omega_{0,:,s}>^2
    %       / (w * Omega_{0,ss}).

    Y = sum(conditional_score.^2, 1) ...
        ./ (w * omega_diag');

    % Avoid log(0) caused by numerical underflow
    Y = max(Y, realmin);

    % Barrier function:
    %
    % f(x) = x - 1 - log(x).

    f_Y = Y - 1 - log(Y);

    % Equation (13) of the paper:
    %
    % T_t =
    % [sum_s f(Y_s) - p*g1(w)]
    % /
    % [g2(w) * sqrt(sum_{s1,s2} R_{s1,s2}^4)].

    statistic = ...
        (sum(f_Y) - p * g1) ...
        / standardization;

    W(u) = statistic;

    % Normally inactive when b = Inf
    if W(u) > b
        break;
    end
end

%% Record the replication index

fid = fopen('tt.txt', 'a');
fprintf(fid, '%g\n', tt);
fclose(fid);

end