function price = priceCap(euribors, resetDates, strikes, vol, dates, discounts)
% priceCap: compute the price of a cap under the Libor Market Model (LMM)
%
% INPUTS:
%   euribors    - vector of forward Euribor rates for each period
%   resetDates  - vector of reset dates
%   strikes     - strike levels of the cap
%   vol         - vector of volatilities corresponding to each strike
%   dates       - dates used in the discount curve
%   discounts   - discount factors corresponding to 'dates'
%
% OUTPUT:
%   price       - total cap price (sum of individual caplets)

% Store today's date as the first reset date
today = resetDates(1);

% Get the number of reset periods
N = length(resetDates)-1;
% Get the number of strikes for which to compute prices
M = length(strikes);

% Initialize price vector with zeros
price = zeros(1,M);

% Loop through each reset period to calculate and sum caplet prices
for i=1:N-1
    % Interpolate the discount factor at the caplet payment date
    B = interpolation(discounts, dates, resetDates(i+2));

    % Price the individual caplet using LMM
    caplet = priceCaplet(euribors(i), strikes, vol, today, resetDates(i+1), resetDates(i+2), B);

    % Add caplet price to total cap price
    price = price + caplet;
end

end