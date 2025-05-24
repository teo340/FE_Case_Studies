function OptionPrice = EuropeanOptionMC(F0,K,B,T,sigma,N,flag)
%European option price with Cox,Ross,Rubinstein tree method
%
%INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number of iterations
% flag:  1 call, -1 put
% OUTPUT
% OptionPrice: price of the option

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

% we generate N samples of a standard normal random variable
% with the command randn(N,1). Then we compute the payoff,
% according to the option type defined above.
prices = payoff(F0*exp(-sigma^2/2*T+sigma*sqrt(T)*randn(N,1)));

% take the mean and discount the quantity
OptionPrice = B*sum(prices)/N;
end %function EuropeanOptionMC