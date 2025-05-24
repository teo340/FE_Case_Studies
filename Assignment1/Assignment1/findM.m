function M = findM(blackPrice, F0, K, B, TTM, sigma, flag, spread)
% This function finds the optimal value of M for the CRR model, iterating
% through every value of M until the error is less than the bid-ask spread.
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% TTM:     time-to-maturity
% sigma: volatility
% flag:  1 call, -1 put
% spread:  bid ask spread
% OUTPUT
% M: optimal value of M

M = 1; %fix M as the lowest value possible assumed by this factor in order not to lose any case
crrPrice = EuropeanOptionCRR(F0,K,B,TTM,sigma,M,flag); 
error = abs(crrPrice-blackPrice);     
while error > spread && M < 2^10
    M = M+1; % increment a counter
    OptionPriceCRR = EuropeanOptionCRR(F0,K,B,TTM,sigma,M,flag); 
    error = abs(OptionPriceCRR-blackPrice);
end 

end     % function findM