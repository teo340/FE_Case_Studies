function OptionPrice = EuropeanOptionMCVR(F0,K,B,T,sigma,N,flag)
% This function computes the price of a European option using 
% a Monte Carlo simulation, where we've applied a variance reduction
% technique called `Antithetic Variables`. This procedure is used to
% reduce the variance of the estimator of the option price.
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number of iterations
% flag:  1 call, -1 put

% depending on the flag, we set a different payoff
% Flag different from 1/-1 are not supported and will
% raise an exception.
if flag == 1
    payoff = @(S) max(S-K,0); % call
elseif flag == -1
    payoff = @(S) max(K-S,0); % put
else
    causeException = MException('MATLAB: wrong flag','Flag is not supported');
    throw(causeException);
end

% We generate two different sample of random variables, where the
% sign of the standard normal random variables has been changed.
f1 = F0*exp(-sigma^2/2*T+sigma*sqrt(T)*randn(N,1));
f2 = F0*exp(-sigma^2/2*T-sigma*sqrt(T)*randn(N,1));
% take the average of the two paths
f = (f1+f2)/2;
prices = payoff(f);

% take the mean and discount the quantity
OptionPrice = B*sum(prices)/N;
end %function EuropeanOptionMC