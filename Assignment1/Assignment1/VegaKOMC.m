function vega = VegaKOMC(F0,K,KO,B,T,sigma,N)
% This function computes the Vega of a Knock-Out option,
% using the Monte Carlo Method.
% INPUTS
% F0: Initial price of the underlying asset
% K: Strike price
% KO: Knock-Out barrier
% B: Discount factor
% T: Maturity
% sigma: Volatility
% N: Number of simulations
% OUTPUT
% vega: Vega of the Knock-Out option

dSigma = 0.01;

% We approximate the derivative of the price with respect to sigma
% using a finite difference scheme. Among the many possibilities, we
% choose the central difference method, because it converges with a quadratic order.
% However, since we all also introducing an additional error approximating the price
% using the MC method, the total order of convergence will be way lower than quadratic.
prevPrice = EuropeanOptionKOMC(F0,K,KO,B,T,sigma-dSigma,N);
nextPrice = EuropeanOptionKOMC(F0,K,KO,B,T,sigma+dSigma,N);

% finite difference: vega = (f(x+dx) - f(x-dx)) / (2*dx) * dx
vega = ((nextPrice - prevPrice) / (2*dSigma)) * dSigma;

end     % function VegaKOMC