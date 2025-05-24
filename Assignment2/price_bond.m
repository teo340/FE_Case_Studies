function bond = price_bond(datesSet, ratesSet)
% INPUT:
%   datesSet : Struct containing maturity dates of deposits and swaps, and
%              settlement & expiry dates of futures
%   ratesSet : Struct containing bid and ask rates of the instruments (e.g., deposits, futures, swaps)
% OUTPUT:
%   bond     : Price of the interbank coupon bond calculated using discounted cash flows

% Set the reference date (settlement date) from datesSet
t0 = datesSet.settlement;

% Calculate the coupon rate as the average of bid and ask rates of 7y swap
coupon = mean(ratesSet.swaps(6,:));

% Bootstraps the discount factor curve from the available market data
[dates, discounts] = bootstrap(datesSet, ratesSet);

% Create a vector of dates that includes the settlement date and the 7-year coupon payment dates
dates = [datetime('2-Feb-2023'); datetime('2-Feb-2024'); dates(13:18)];

% Select the discount factors from the bootstrap that correspond to the 7-year coupon payments
discounts = [0.968301310232409; discounts(13:18)];

% Compute the year fraction between the settlement date and each cash
% flow date using day count convention 6 (30/360 European)
yf = yearfrac(t0, dates, 6);

% Calculate the time intervals (in years) between successive cash flow dates
deltas = yf(2:end)-yf(1:end-1);

% Calculate the bond price:
% - Each coupon payment is discounted
% - The bond price is the sum of these discounted coupons plus the discounted nominal value (1) at maturity
bond = coupon*sum(deltas.*discounts) + discounts(end);
end