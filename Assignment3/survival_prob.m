function p = survival_prob(T)
% This function computes the survival probability at time T.
%
% This function models a 'piecewise' exponential survival probability,
% where the default intensity changes at a threshold time theta.
%
% Inputs:
% T: the time at which the survival probability is evaluated.
%
% Outputs:
% p: the survival probability at time T.

% Define Parameters
bp = 1e-4;           % Basis point
lambdas = [5;9]*bp;  % Intensities of the model (before and after theta)
theta = 4;           % Time threshold where intensity changes

% Compute survival probability
if T<theta
    % Before threshold theta, survival follows exp(-lambda1 * T)
    p = exp(-lambdas(1)*T);
else
    % After theta, survival follows exp(-lambda1*theta - lambda2*(T-theta))
    p = exp(-lambdas(1)*theta-lambdas(2)*(T-theta));
end

end
