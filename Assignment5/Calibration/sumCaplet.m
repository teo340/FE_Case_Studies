function price = sumCaplet(euribors, strike, vol, resetDates, today, discounts)
% sumCaplet - Calculates the sum of individual caplets
%
% Inputs:
%   euribors - Vector of forward Euribor rates
%   strike - Strike rate for the cap
%   vol - Vector of volatilities corresponding to each caplet
%   resetDates - Vector of reset dates for caplets
%   today - Valuation date

%   discounts - Vector of discount factors for payment dates
%
% Output:
%   price - Sum of all caplet prices

% Get the number of caplets to price
N = length(euribors);

% Initialize price to zero
price = 0;

% Loop through each caplet and sum their prices
for i = 1:N
    price = price + priceCaplet(euribors(i), strike, vol(i,:), today, resetDates(i), resetDates(i+1), discounts(i));
end

end