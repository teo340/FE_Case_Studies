function optionPrice=EuropeanOptionKOMC(F0,K,KO,B,T,sigma,N)
% This function calculates the price of a European option with a knock-out,
% using a Monte Carlo approach. 
% INPUTS
% F0: current price of the underlying asset
% K: strike price
% KO: knock-out barrier
% B: discount factor
% T: time to maturity
% sigma: volatility of the underlying asset
% N: number of simulations
% OUTPUT
% optionPrice: price of the European option with knock-out

% Since the option is European, it is enough to check if the
% barrier has been hit at maturity. If so, the payoff is set to zero.
% Otherwise, the payoff is the same as for a standard European option.
% We simulate N paths of the underlying asset price and peform the
% barrier check for each path.
prices = zeros(N,1);
for i=1:N
    g = randn();
    s = F0*exp(-sigma^2/2*T+sigma*sqrt(T)*g);
    if s > KO
        % if the barrier has been hit, the payoff is zero
        prices(i) = 0;
    else
        prices(i) = max(F0*exp(-sigma^2/2*T+sigma*sqrt(T)*g)-K,0);
    end
end
    
% Sample mean of the payoffs, discounted
optionPrice = B*sum(prices)/N;
end