function optionPrice = EuropeanOptionKOCRR(F0,K,KO,B,T,sigma,N)
% We compute the price of the Knock-Out Call option using
% a tree-based approach, based on the Cox-Ross-Rubinstein model.
% INPUTS
% F0: initial stock price
% K: strike price
% KO: knock-out barrier
% B: discount factor
% T: maturity
% sigma: volatility
% N: number of time steps
% OUTPUT
% optionPrice: price of the European option with knock-out

dt = T/N; % time step
u = exp(sigma*sqrt(dt)); % up factor
d = 1/u; % down factor
q = (1-d)/(u-d); % probability of up movement

r = -log(B) / T;
% discount factor for the time interval
B_dt = exp(-r*dt);

% compute leaves of the tree
leaves = F0*u.^(-N:2:N);
payoff = zeros(N+1,1);

% For a European Knock-Out Call option, it is enough to check
% if the barrier has been hit at maturatiy. If so,
% the payoff is set to zero. Then we procede with the standard
% procedure to price European options with a Tree
for i=1:length(leaves)
    if leaves(i) > KO
        % if the barrier is hit, the payoff is zero
        payoff(i) = 0;
    else
        % otherwise, classical payoff
        payoff(i) = max(leaves(i)-K,0);
    end
end

% Standard approach for trees
for i=N:-1:1
    for j=1:i
        payoff(j) = B_dt*((1-q)*payoff(j)+q*payoff(j+1));
    end
end

optionPrice = B_dt*payoff(1);
end     % function EuropeanOptionKOCRR