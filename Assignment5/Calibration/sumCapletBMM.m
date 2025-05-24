function price = sumCapletBMM(euribors, strike, v, resetDates, today, discounts)
% sumCapletBMM - Calculates the total price of a cap as the sum of individual caplets
% using the Bond Market Model (BMM) approach
%
% Inputs:
%   euribors - Vector of forward Euribor rates
%   strike - Strike rate for the cap
%   v - Vector of BMM volatilities corresponding to each caplet
%   resetDates - Vector of reset dates for caplets
%   today - Valuation date
%   discounts - Vector of discount factors for payment dates
%
% Output:
%   price - Total price of the cap (sum of all caplet prices under BMM)

% Get the number of reset dates
N = length(resetDates);

% Initialize price to zero
price = 0;

% Loop through each caplet period (number of caplets = N-1)
for i = 1:N-1
    price = price + priceCapletBMM(euribors(i), strike, v(i,:), today, resetDates(i), resetDates(i+1), discounts(i));
end

end