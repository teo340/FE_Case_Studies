function vega = VegaKOCRR(F0,K,KO,B,T,sigma,N)
% This function computes the Vega of a Knock-Out option,
% using the CRR approach.
% INPUTS
% F0: Initial price of the underlying asset
% K: Strike price
% KO: Knock-Out barrier
% B: Discount factor
% T: Maturity
% sigma: Volatility
% N: Number of steps for CRR
% OUTPUT
% vega: Vega of the Knock-Out option

dSigma = 0.01; % 1% of basis point

% We approximate the derivative of the price with respect to sigma
% using a finite difference scheme. Among the many possibilities, we
% choose the central difference method, because it converges with a quadratic order.
% However, since we all also introducing an additional error approximating the price
% using the CRR method, the total order of convergence will be way lower than quadratic.
prevPrice = EuropeanOptionKOCRR(F0,K,KO,B,T,sigma-dSigma,N);
nextPrice = EuropeanOptionKOCRR(F0,K,KO,B,T,sigma+dSigma,N);

% finite difference: vega = (f(x+dx) - f(x-dx)) / (2*dx) * dx
vega = ((nextPrice - prevPrice) / 2);

end     % function VegaKOCRR