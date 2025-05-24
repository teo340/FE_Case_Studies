function [lambdas_est, CI] = estimate_params(tau_samples, theta)
% This function estimates the intensity parameters lambda1 and lambda2 
% from simulated default times and computes their confidence intervals.
%
% Inputs:
% tau_samples: a vector of simulated default times.
% theta: the threshold time at which the default intensity changes.
%
% Outputs:
% lambdas_est: estimated values of lambda1 and lambda2.
% CI: confidence intervals for lambda1 and lambda2.
M = length(tau_samples); % Number of samples
counters = zeros(30,1);
% Count how many samples have survived past each time point
for i=1:30
    for j=1:M
        if tau_samples(j) >= i
            counters(i) = counters(i)+1;
        end
    end
end

% Compute simulated survival probabilities
prop = counters/M;
% Compute theoretical survival probabilities using the function 'survival_prob'
survs = zeros(30,1);
for i=1:30
    survs(i) = survival_prob(i);
end

% Plot theoretical and simulated survival probabilities
plot(1:30, survs, '--');
grid on; hold on;
plot(1:30, prop(1:30));
legend('Theoretical','Simulated')

% plot(1:30, survs, '--');
% grid on; hold on;
% plot(1:30, prop(1:30));

% Now we can proceed with the separation the data and with the estimation of the lambda parameters
% Fit the model to the first 4 data points (lambda1)
log_p1 = -log(prop(1:4));             % Log transformation
lam1_est = (1:4)'\log_p1;             % Estimate lambda1 using least squares

% Fit the model to data from time 5 to 30 (lambda2)
log_p2 = log(prop(5:30));              % Log transformation
X = [ones(length(log_p2),1),(5:30)'];  % Regression matrix
beta = X\log_p2; % Linear regression

% Store estimated parameters
lambdas_est(1) = lam1_est;
lambdas_est(2) = -beta(2); 

% From now on we proceed with the computation of the confidence Intervals
y = log_p2;
y_hat = X*beta;
residuals = y - y_hat;     % Residuals
n = size(X, 1);            % Number of observations
p = size(X, 2);            % Number of regressors (including intercept if present)

sigma2 = sum(residuals.^2) / (n - p); % Estimate of error variance

C = inv(X' * X);                   % Covariance matrix of beta
se_beta = sqrt(diag(sigma2 * C));  % Standard errors

% Compute 95% confidence intervals using t-distribution
alpha = 0.05; 
t_crit = tinv(1 - alpha/2, n - p);  % Critical t-value

% Compute confidence intervals for estimated parameters
CI_lower = beta - t_crit * se_beta;
CI_upper = beta + t_crit * se_beta;

% Compute confidence intervals for lambda1 and lambda2
CI_lambda1 = lambdas_est(2)-[CI_upper(1),CI_lower(1)]/theta;
CI_lambda2 = [-CI_upper(2),-CI_lower(2)];

% Store confidence intervals
CI = [CI_lambda1; CI_lambda2];   

end
