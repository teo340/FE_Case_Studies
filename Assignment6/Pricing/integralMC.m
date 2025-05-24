function [price, IC] = integralMC(sigma, k, eta, t, logMoneyness, F0, discount_1y)
% This function numerically integrates the Lewis Formula using Monte Carlo
% INPUTS:
% eta: skewness
% sigma: level of volatility
% k: volatility of volatility
% t: maturity date
% logMoneyness: vector of logMoneyness values
% F0: forward price
% discount_1y: discount factor
% OUTPUTS:
% price: vector of call option prices for each logMoneyness value
% IC: confidence interval for the price

% Monte Carlo parameters
N_sim = 1e5;
rng(42);
g = randn(N_sim,1);

a = t/k; b = 1/a;

G = gamrnd(a,b,[N_sim,1]);

% To assess to goodness of the generate Gamma samples, 
% we compute the first four moments of the Gamma distribution
% and compare them with the empirical moments of the generated samples.
mu1 = mean(G);
mu2 = mean(G.^2);
mu3 = mean(G.^3);
mu4 = mean(G.^4);

fprintf("1st Moment: %0.5f - %d\n",mu1,a*b);
fprintf("2nd Moment: %0.5f - %d\n",mu2,a*b^2*(a+1));
fprintf("3rd Moment: %0.5f - %d\n",mu3,a*b^3*(a+1)*(a+2));
fprintf("4th Moment: %0.5f - %d\n",mu4,a*b^4*(a+1)*(a+2)*(a+3));

% Characteristic function of the forward price process
laplaceExponent = @(u) -t/k*log(1+k*u*sigma^2);
mu = -laplaceExponent(eta);
f = mu - (0.5+eta)*sigma^2*G*t+sigma*sqrt(t)*sqrt(G).*g;
Ft = F0*exp(f);
strike = F0./exp(logMoneyness);

payoff = discount_1y * max(Ft - strike,0);
price = mean(payoff);

% Confidence interval
alpha = 0.95;
stdev = std(payoff);

IC = [price - norminv(1-alpha/2)*stdev; price + norminv(1-alpha/2)*stdev];

end