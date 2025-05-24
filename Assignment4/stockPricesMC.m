function [S1,S2] = stockPricesMC(rho, sigma1, sigma2, d1, d2, S0, startDate, discounts, dates, M)
% This functions simulates the paths of two correlated stocks,
% according to a Black dynamics. Prices are computed in each ref date.
% INPUTS:
% rho: correlation
% sigma1,sigma2: volatility
% d1,d2: dividend yield
% S0: initial prices, vector
% startDate: start date
% discounts,dates: from the boostrap, used for interpolation
% M: number of simulations
% OUTPUTS:
% S1,S2: simulated prices of the two stocks


ACT_365 = 3;

% reference dates
ref_dates = businessdayoffset(startDate+calyears(0:4));
deltas = yearfrac(ref_dates(1:end-1), ref_dates(2:end), ACT_365);

% interpolated discounts and fwd discounts
ref_discounts = [1 interpolation(discounts, dates, ref_dates(2:end))];
fwd_discounts = ref_discounts./[1 ref_discounts(1:end-1)];

S1 = ones(M,length(ref_dates)); S1(:,1) = S0(1);
S2 = ones(M,length(ref_dates)); S2(:,1) = S0(2);

for i=1:length(ref_dates)-1
    % Multivariate gaussian random variables
    Z = mvnrnd([0,0], [1,rho; rho,1], M);
    Z1 = Z(:,1); Z2 = Z(:,2);
    
    S1(:,i+1) = S1(:,i)/fwd_discounts(i+1).*exp((-d1-0.5*sigma1^2)*deltas(i)+sigma1*sqrt(deltas(i))*Z1);
    S2(:,i+1) = S2(:,i)/fwd_discounts(i+1).*exp((-d2-0.5*sigma2^2)*deltas(i)+sigma2*sqrt(deltas(i))*Z2);
end

end