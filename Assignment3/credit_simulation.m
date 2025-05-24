function [tau_samples] = credit_simulation(lambdas, theta, M)
% This function simulates default times based on a two-phase intensity model.
%  lambdas: a vector containing two intensity parameters [lambda1, lambda2].
%  theta: a threshold time at which the intensity changes.
%  M: the number of default times to simulate.
%
%  Output:
%  tau_samples: a vector (of dimension M) containing the simulated default times.

% Set the random number generator to its default state for reproducibility
rng('default');
% Generate M uniform random numbers between 0 and 1
u = rand(M,1);
% Initialize the output vector with zeros 
tau_samples = zeros(M,1);
% Loop through M simulations
for i=1:M
    % Compute tau according to the theory
    tau = -log(1-u(i))/lambdas(1);
    % Check if the generated time is below the threshold
    if tau < theta
     % If yes, assign the generated value directly   
        tau_samples(i) = tau;
    else
     % If tau exceeds theta, adjust it using the second intensity lambda2
        tau_samples(i) = (-log(1-u(i))+lambdas(2)*theta-lambdas(1)*theta)/lambdas(2);
    end
end
end
