function optionPrice = EuropeanOptionCRR(F0,K,B,T,sigma,N,flag)
%European option price with Cox,Ross,Rubinstein tree method
%
%INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number of time intervals
% flag:  1 call, -1 put

dt = T/N; % time step
u = exp(sigma*sqrt(dt)); % up factor
d = 1/u; % down factor
q = (1-d)/(u-d); % probability of up movement

r = -log(B) / T;
% discount factor for the time interval
B_dt = exp(-r*dt);

% Payoff can either be Call or Put.
% Other option types are not supported.

if flag == 1 % call option
    payoff = max(F0*u.^(-N:2:N)-K,0);
elseif flag == -1 % put option
    payoff = max(K-F0*u.^(-N:2:N),0);
else
    causeException = MException('MATLAB: wrong flag','Flag is not supported');
    throw(causeException);
end

% We compute the payoff through a backward iteration,
% from the leaves of the tree to the root.

for i = N:-1:1
    % vectorization is faster
    payoff(1:i) = B_dt * ((1 - q) * payoff(1:i) + q * payoff(2:i+1));
end

optionPrice = B_dt*payoff(1);

end %function EuropeanOptionCRR