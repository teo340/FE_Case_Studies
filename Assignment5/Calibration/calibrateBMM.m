function [vBMM] = calibrateBMM(euribors, strikes, resetDates, vol, dates, discounts)
% calibrateBMM - Calibrates a Bond Market Model given market flat volatilities
%
% Inputs:
%   euribors - Forward Euribor rates
%   strikes - Strike prices for caps
%   resetDates - Reset dates for caplets
%   vol - Initial volatility matrix (by tenor and strike)
%   dates - Dates for discount factors
%   discounts - Discount factors
%
% Output:
%   vBMM - Calibrated Bond Market Model v-matrix (by reset date and strike)

today = dates(1);

% Get dimensions of the vol matrix
[~, m] = size(vol);

% Initialize BMM v-matrix for 39 reset dates (4*10-1) and m strikes
vBMM = zeros(4*10-1, m);

ACT_365 = 3;

capBMM = cell(1, length(strikes));

% Calibrate first reset date (3-month caplet)
cap = priceCap(euribors(1), resetDates(1:3), strikes, vol(1,:), dates, discounts);
for i = 1:length(strikes)
    % Get discount factors for this caplet
    discountsCap = interpolation(discounts,dates, resetDates(3));
    
    % Create function to calculate BMM caplet price with volatility parameter v
    capBMM{i} = @(v) sumCapletBMM(euribors(1), strikes(i), v, resetDates(2:3), resetDates(1), discountsCap);

    % Create objective function: difference between BMM price and market price
    minimizer = @(v) capBMM{i}(v) - cap(i);

    % Find volatility that matches market price using root-finding
    vBMM(1, i) = fzero(minimizer, 0.1);
end

% Calibrate second reset date (6-month caplet)
cap = priceCap(euribors(1:2), resetDates(1:4), strikes, vol(1,:), dates, discounts);
for i = 1:length(strikes)
    % Get discount factors for this caplet
    discountsCap = interpolation(discounts,dates, resetDates(3:4));
    
    % Create function to calculate BMM caplet price with volatility parameter v
    % Note: Using previously calibrated volatility for first reset date
    capBMM{i} = @(v) sumCapletBMM(euribors(1:2), strikes(i), [vBMM(1,i); v], resetDates(2:4), resetDates(1), discountsCap);

    % Create objective function: difference between BMM price and market price
    minimizer = @(v) capBMM{i}(v) - cap(i);

    % Find volatility that matches market price using root-finding
    vBMM(2, i) = fzero(minimizer, 0.1);
end

% Calibrate third reset date (9-month caplet)
cap = priceCap(euribors(1:3), resetDates(1:5), strikes, vol(1,:), dates, discounts);
for i = 1:length(strikes)
    % Get discount factors for this caplet
    discountsCap = interpolation(discounts,dates, resetDates(3:5));
    
    % Create function to calculate BMM caplet price with volatility parameter v
    % Note: Using previously calibrated volatilities for first two reset dates
    capBMM{i} = @(v) sumCapletBMM(euribors(1:3), strikes(i), [vBMM(1,i); vBMM(2,i); v], resetDates(2:5), resetDates(1), discountsCap);

    % Create objective function: difference between BMM price and market price
    minimizer = @(v) capBMM{i}(v) - cap(i);

    % Find volatility that matches market price using root-finding
    vBMM(3, i) = fzero(minimizer, 0.1);
end

% Store 1-year cap prices for next iteration
capPrev = cap;

% Loop through remaining tenors (2 to 10 years)
for j = 2:10
    % Price cap up to current tenor using LMM for reference
    capNext = priceCap(euribors(1:4*j-1), resetDates(1:4*j+1), strikes, vol(j,:), dates, discounts);
    
    % Calculate incremental cap price difference
    deltaC = capNext - capPrev;
    
    % Interpolate discount factors for the relevant reset dates
    discountsCap = interpolation(discounts, dates, resetDates(4*j-2:4*j+1));
    
    % Calibrate for each strike
    for i = 1:m
        % Create linear interpolation functions for intermediate v
        % First intermediate point
        vST = @(vTH) vBMM(4*(j-1)-1, i) + (vTH - vBMM(4*(j-1)-1, i)) / ...
              yearfrac(resetDates(4*(j-1)+1), resetDates(4*j+1), ACT_365) * ...
              yearfrac(resetDates(4*(j-1)+1), resetDates(4*(j-1)+2), ACT_365);
        
        % Second intermediate point
        vND = @(vTH) vBMM(4*(j-1)-1, i) + (vTH - vBMM(4*(j-1)-1, i)) / ...
              yearfrac(resetDates(4*(j-1)+1), resetDates(4*j+1), ACT_365) * ...
              yearfrac(resetDates(4*(j-1)+1), resetDates(4*(j-1)+3), ACT_365);
        
        % Third intermediate point
        vRD = @(vTH) vBMM(4*(j-1)-1, i) + (vTH - vBMM(4*(j-1)-1, i)) / ...
              yearfrac(resetDates(4*(j-1)+1), resetDates(4*j+1), ACT_365) * ...
              yearfrac(resetDates(4*(j-1)+1), resetDates(4*(j-1)+4), ACT_365);
        
        % Create objective function: difference between BMM price and market price difference
        f = @(vTH) sumCapletBMM(euribors(4*(j-1):4*j-1), strikes(i), ...
            [vST(vTH); vND(vTH); vRD(vTH); vTH], ...
            resetDates(4*(j-1)+1:4*j+1), today, discountsCap) - deltaC(i);
        
        % Find volatility that matches market price difference using root-finding
        vBMM(4*j-1, i) = fzero(f, vBMM(4*j-5, i));
        
        % Assign interpolated v to intermediate reset dates
        vBMM(4*j-2, i) = vRD(vBMM(4*j-1, i));
        vBMM(4*j-3, i) = vND(vBMM(4*j-1, i));
        vBMM(4*j-4, i) = vST(vBMM(4*j-1, i));
    end
    
    % Update previous cap price for next iteration
    capPrev = capNext;
end

end