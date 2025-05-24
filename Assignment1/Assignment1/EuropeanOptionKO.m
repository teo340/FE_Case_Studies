function optionPrice = EuropeanOptionKO(F0,K,KO,B,T,sigma)
% This function calculates the price of a European option with a knock-out,
% with a closed formula. The Knock-out barrier can be replicated with a
% portfolio of long Bull spread and short (K-KO) digital option.
% INPUTS
% F0: current price of the underlying asset
% K: strike price
% KO: knock-out barrier
% B: discount factor
% T: time to maturity
% sigma: volatility of the underlying asset
% OUTPUT
% optionPrice: price of the European option with knock-out

% We compute the price of the bull spread (CallK - CallKO)
% using the previoulsy defined function EuropeanOptionClosed
callK = EuropeanOptionClosed(F0,K,B,T,sigma,1);
callKO = EuropeanOptionClosed(F0,KO,B,T,sigma,1);

% We compute the price of the digital option using a closed formula
% derived through direct calculation, according to Black Model
d2 = log(F0/KO)/(sigma*sqrt(T))-0.5*sigma*sqrt(T);
digital = (KO-K)*B*normcdf(d2);

% The sum of the Knock-out price is the sum of the components
optionPrice = callK - callKO - digital;

end     %function EuropeanOptionKO