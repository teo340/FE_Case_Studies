function [M,errorCRR]=PlotErrorCRR(F0,K,B,T,sigma)
% Plot the errors of the Cox, Ross, Rubinstein method
% for the pricing of a European Option.
%
% INPUT:
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     expiry
% sigma: volatility
% OUTPUT:
% M:     number of steps
% errorCRR: error of the CRR method

% benchmark, price computed with the closed formula
actualPrice = EuropeanOptionClosed(F0,K,B,T,sigma,1);

% we set up several binomial trees, each with M(i) = 2^i
% time steps. 
m=1:10;
M=2.^m;

% For each tree, we compute the error, in absolute value,
% with respect to the price obtained with the closed formula.
errorCRR = zeros(1,length(M));
for i=1:length(M)
    crrPrice = EuropeanOptionCRR(F0,K,B,T,sigma,M(i),1);
    errorCRR(i) = abs(actualPrice-crrPrice);
end

spread = ones(1,length(M))*1e-4; % spread 1bp, target error

figure;
loglog(M,errorCRR,'LineWidth',1.5);
hold on; grid on;
loglog(M,1./M,'LineWidth',1.5);
loglog(M,spread,'LineWidth',1.5,'LineStyle','--');
title('Error CRR');
legend('abserr','1/M','spread');
hold off;



end     %function PlotErrorCRR