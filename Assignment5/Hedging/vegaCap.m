function vega = vegaCap(euribors, mkt_strikes, resetDates, sigma, dates, discounts, swapRate)
% vegaCap computes the total vega of a cap by summing the vega of its individual caplets.
%
% INPUTS:
%   euribors     - Vector of Euribors
%   mkt_strikes  - Vector of market strikes used for volatility interpolation
%   resetDates   - Vector of reset dates 
%   sigma        - Matrix of volatilities 
%   dates        - Dates corresponding to the discount curve
%   discounts    - Discount factors corresponding to 'dates'
%   swapRate     - Strike rate of the cap
%
% OUTPUT:
%   vega         - Total Black vega of the cap (sensitivity to volatility)

% Day count conventions
ACT_360 = 2;  
ACT_365 = 3;  

% Initialize vega 
vega = 0;

% Loop over each caplet (note: ends at length-2 to access i+2 safely)
for i = 1:length(resetDates) - 2
    % Interpolate discount factor 
    B = interpolation(discounts, dates, resetDates(i+2));

    % Year fraction 
    yf = yearfrac(resetDates(i+1), resetDates(i+2), ACT_360);

    % Time to expiry from today to the start of the caplet
    deltaT = yearfrac(resetDates(1), resetDates(i+1), ACT_365);

    % Interpolate volatility for the given swapRate from the market vol surface
    int_vol = spline(mkt_strikes, sigma(i,:), swapRate);

    % Compute d1 term in the Black model for this caplet
    d1 = log(euribors(i) / swapRate) / (int_vol * sqrt(deltaT)) + 0.5 * int_vol * sqrt(deltaT);

    % Accumulate the Black vega of this caplet
    vega = vega + euribors(i) * yf * B * sqrt(deltaT) * exp(-d1^2 / 2) / sqrt(2 * pi);
end

end
