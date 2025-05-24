function [upfront, CI] = computeUpfront3y(spread, dates, discounts)
% The function computes the upfront payment for a 3-year product
% using Monte Carlo simulation. Differently fromt he 2-year product,
% the path dependency of the payoff makes it impossible to price the certificate
% with a closed formula.
% INPUTS:
% spread: spread of the certificate
% dates: vector of dates from boostrap, datetime
% discounts: vector of discount factors from bootstrap

% Variance-Gamma parameters
sigma = 0.2;
k = 1;
eta = 3;
alpha = 0; % Variance Gamma

% Day count conventions
ACT_360 = 2;
ACT_365 = 3;

today = datetime("02-Feb-2023");
maturity_1y = businessdayoffset(today + calyears(1));
payment_1y = businessdayoffset(maturity_1y + caldays(2));
maturity_2y = businessdayoffset(today + calyears(2));
payment_2y = businessdayoffset(maturity_2y + caldays(2)); 
maturity_3y = businessdayoffset(today + calyears(3));
payment_3y = businessdayoffset(maturity_3y + caldays(2)); 

% Relevant market data
load("./Dataset/cSelect.mat");
S0 = cSelect.reference;
d = cSelect.dividend;
T_1y = yearfrac(today, maturity_1y, ACT_365);

B_1y = interpolation(discounts, dates, payment_1y);
B_2y = interpolation(discounts, dates, payment_2y);
B_3y = interpolation(discounts, dates, payment_3y);

% Compute the forward price F0
F0_1y = S0 / B_1y * exp(-d * T_1y);

strike = 3200;

resetDates = businessdayoffset(today:calmonths(3):maturity_3y);
paymentDates = businessdayoffset(resetDates+caldays(2));

yf = yearfrac(paymentDates(1:end-1), paymentDates(2:end), ACT_360);
discounts_inter = interpolation(discounts, dates, paymentDates(2:end));
BPV = sum(yf.*discounts_inter);

% The coupon is simulated using the Monte Carlo method.
[coupon,Ft] = couponMC(F0_1y, strike, sigma, k, eta, today, maturity_3y, [B_1y, B_2y, B_3y]);

BPV_1y = sum(yf(1:4).*discounts_inter(1:4));
BPV_2y = sum(yf(5:8).*discounts_inter(5:8));
BPV_3y = sum(yf(9:end).*discounts_inter(9:end));

prob_up = Ft(:,2) >= strike;
prob_up2 = Ft(:,2) >= strike & Ft(:,3) >= strike;

upfront = 1 - B_1y + spread*BPV_1y + ...
          prob_up*(spread*BPV_2y - B_2y + B_1y) + ...
          prob_up2*(spread*BPV_3y - B_3y + B_2y) - coupon;

N_sim = length(upfront);
sd = std(upfront);
upfront = mean(upfront);

% we provide a confidence interval alpha = 0.05 for the upfront
CI = [upfront - norminv(0.975)*sd/N_sim, upfront + norminv(0.975)*sd/N_sim];

end