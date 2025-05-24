function [M,stdEstim]=PlotErrorMC(F0,K,B,T,sigma)
% Plot the errors of the MonteCarlo method
% for the pricing of an European Option.
%
% INPUT:
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     expiry
% sigma: volatility
% OUTPUT:
% M:         number of simulations
% stdEstim:  standard deviation of the estimations

% we set up 20 different montecarlo routines,
% each with 2^i number of iterations
m = 1:20;
M = 2.^m;
stdEstim = zeros(length(M),1);

% we generate M(i) samples of a standard gaussian random variable, 
% to then compute the discounted payoff.
for i=1:length(M)
    prices = B*max(F0*exp(-sigma^2/2*T+sigma*sqrt(T)*randn(M(i),1))-K,0);
    stdEstim(i) = 1/sqrt(M(i))*std(prices);
end

spread = ones(1,length(M))*1e-4; % spread is 1bd, target for our method

figure;
loglog(M,stdEstim,'LineWidth',1.5);
hold on; grid on;
loglog(M,1./sqrt(M),'LineWidth',1.5);
loglog(M,spread,'LineWidth',1.5,'LineStyle','--');
title('Error MC');
legend('stddev','1/sqrt(M)','spread');
hold off;

end     %function PlotErrorMC