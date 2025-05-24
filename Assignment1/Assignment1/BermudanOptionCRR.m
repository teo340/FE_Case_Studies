function optionPriceBerm = BermudanOptionCRR(F0, K, B, T, sigma, d, N)
% Bermudan option price with CRR method and Pseudo-American option price
%
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% d:     dividend yield
% N:     number of steps
%
% OUTPUT
% optionPriceBerm: price of the Bermudan option with 2 exercise dates

% reduce N to be a multiple of 3
N = 4*floor(N/4);
% dt is the interval of time between two knots
dt = T/N;
% u is the factor by which the stock price increases
u = exp(sigma * sqrt(dt));
% q is the probability under the risk neutral measure
q = 1 / (u + 1);

% Compute the discount factor for each step
r = -log(B) / T;
B_dt = exp(-r*dt);

% initialize the tree leaves (N+1) with the forward and payoff
Ftt = F0 * u.^(-N:2:N);
leavesCRR = max(Ftt - K, 0);

% reduce the tree to the root for the Bermudan option
for i = N-1:-1:0
    % update the option price
    leavesCRR = B_dt * ((1-q)*leavesCRR(1:end-1) + q*leavesCRR(2:end));
    if i == N/4 || i == N/2 || i == 3*N/4
        % compute the price of the forward at time i
        Fti = F0 * u.^(-i:2:i);
        % We discount using also the dividends
        Sti = Fti * exp(-((r-d)*(N-i)*dt));
        % Instead of comparing the max twice, we can just compare once
        leavesCRR = max(leavesCRR, Sti - K); 
    end
end

optionPriceBerm = leavesCRR;
end     % function BermudanOptionCRR