function [upfront] = computeUpfront(spread, dates, discounts, flag)
% This function computes the upfront payment of a certificate for a given spread.
% The flag parameter is used to control the pricing model used:
%  - "VG" for Variance Gamma
%  - "Black" for Black model
%  - "Smile" for Black with smile adjustment for digital risk
% INPUTS:
% spread: spread of the certificate
% dates: vector of dates from boostrap, datetime
% discounts: vector of discount factors from bootstrap
% flag: pricing model to use ('VG', 'Black', 'Smile')

% Variance-Gamma parameters
sigma = 0.2;
k = 1;
eta = 3;
alpha = 0;

% Day count conventions
ACT_360 = 2;
ACT_365 = 3;
EU_30_360 = 6;

today = dates(1);
maturity_1y = businessdayoffset(today + calyears(1));
payment_1y = businessdayoffset(maturity_1y + caldays(2));
maturity_2y = businessdayoffset(today + calyears(2));
payment_2y = businessdayoffset(maturity_2y + caldays(2)); 

% Relevant market data
load("./Dataset/cSelect.mat");
mkt_strikes = cSelect.strikes;
mkt_vols = cSelect.surface;
S0 = cSelect.reference;
d = cSelect.dividend;
T_1y = yearfrac(today, maturity_1y, ACT_365);

B_1y = interpolation(discounts, dates, payment_1y);
B_2y = interpolation(discounts, dates, payment_2y);

% Compute the forward price F0
F0_1y = S0 / B_1y * exp(-d * T_1y);

strike = 3200;

resetDates = businessdayoffset(today:calmonths(3):maturity_2y);
paymentDates = businessdayoffset(resetDates+caldays(2));

yf = yearfrac(paymentDates(1:end-1), paymentDates(2:end), ACT_360);
discounts_inter = interpolation(discounts, dates, paymentDates(2:end));

BPV_1y = sum(yf(1:4).*discounts_inter(1:4));
BPV_2y = sum(yf(5:end).*discounts_inter(5:end));

delta_1y = yearfrac(today, maturity_1y, EU_30_360);
delta_2y = yearfrac(maturity_1y, maturity_2y, EU_30_360);

% The upfront is determined by computing the values of the two coupons.
% Note that, while the first coupon will always be paid, the second coupon
% will be paid only if the early redemption clause is not triggered.
% Hence, we compute the NPV of the two coupons and the probability
% prob_up of not triggering the early redemption clause.

if flag == "VG"
    prob_up = digitalPrice(strike, sigma, k, eta, alpha, F0_1y, 1, T_1y);

    coupon_1y = 0.06*delta_1y*B_1y*(1-digitalPrice(strike, sigma, k, eta, alpha, F0_1y, 1, T_1y));
    coupon_2y = 0.02*delta_2y*B_2y*digitalPrice(strike, sigma, k, eta, alpha, F0_1y, 1, T_1y);

    upfront = 1 - B_1y + spread*BPV_1y + prob_up*(spread*BPV_2y - B_2y + B_1y) - (coupon_1y+coupon_2y);
elseif flag == "Black"
    prob_up = digitalPriceBlk(mkt_strikes, mkt_vols, strike, F0_1y, 1, T_1y, 'BLK');

    coupon_1y = 0.06*delta_1y*B_1y*(1-digitalPriceBlk(mkt_strikes, mkt_vols, strike, F0_1y, 1, T_1y, 'BLK'));
    coupon_2y = 0.02*delta_2y*B_2y*digitalPriceBlk(mkt_strikes, mkt_vols, strike, F0_1y, 1, T_1y, 'BLK');

    upfront = 1 - B_1y + spread*BPV_1y + prob_up*(spread*BPV_2y - B_2y + B_1y) - (coupon_1y+coupon_2y);
elseif flag == "Smile"
    prob_up = digitalPriceBlk(mkt_strikes, mkt_vols, strike, F0_1y, 1, T_1y, 'Smile');

    coupon_1y = 0.06*delta_1y*B_1y*(1-digitalPriceBlk(mkt_strikes, mkt_vols, strike, F0_1y, 1, T_1y, 'Smile'));
    coupon_2y = 0.02*delta_1y*B_2y*digitalPriceBlk(mkt_strikes, mkt_vols, strike, F0_1y, 1, T_1y, 'Smile');

    upfront = 1 - B_1y + spread*BPV_1y + prob_up*(spread*BPV_2y - B_2y + B_1y) - (coupon_1y+coupon_2y);
else
    error("Flag not defined")
end
