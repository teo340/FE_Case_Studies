function [sigma] = calibrateLMM_full(euribors, strikes, resetDates, vol, dates, discounts)
% This function calibrates the LMM model to the given market data, using every date.
% It assumes that data is yearly up to 10 years, afterwards it follows the same convention
% of the table provided (12y, 15y, 20y).
% INPUTS:
% euribors: 3m forward libor rates
% strikes
% resetDates: reset dates of the caplets, datetime
% vol: matrix of volatilities to calibrate on
% dates: bootstrap dates, datetime
% discounts: bootstrap discounts

today = dates(1);
[~, m] = size(vol);
sigma = zeros(79, m);

sigma(1:3, :) = repmat(vol(1, :), 3, 1);

cap1y = priceCap(euribors(1:3), resetDates(1:5), strikes, vol(1, :), dates, discounts);
capPrev = cap1y;

% Calibration in the first 10 years
for j = 2:10
    capNext = priceCap(euribors(1:4*j-1), resetDates(1:4*j+1), strikes, vol(j, :), dates, discounts);
    deltaC = capNext - capPrev;
    discountsCap = interpolation(discounts, dates, resetDates(4*j-2:4*j+1));
    
    for i = 1:m
        vec = @(sigmaTH) sigmaVec(sigmaTH, sigma, resetDates, j, i, 4);
        % function to minimize
        f = @(sigmaTH) sumCaplet(euribors(4*(j-1):4*j-1), strikes(i), vec(sigmaTH), resetDates(4*(j-1)+1:4*j+1), today, discountsCap) - deltaC(i);
        sigma(4*j-1, i) = fzero(f, 0.2);
        tmp = vec(sigma(4*j-1, i));
        for k = 1:3
            sigma(4*j-k-1, i) = tmp(end-k);
        end
    end
    
    capPrev = capNext;
end

% Calibration from 10 to 12 years
capNext = priceCap(euribors(1:47), resetDates(1:49), strikes, vol(11, :), dates, discounts);
deltaC = capNext - capPrev;
discountsCap = interpolation(discounts, dates, resetDates(41:49));

for i = 1:m
    vec = @(sigmaTH) sigmaVec(sigmaTH, sigma, resetDates, 11, i, 8);
    % function to minimize
    f = @(sigmaTH) sumCaplet(euribors(40:47), strikes(i), vec(sigmaTH), resetDates(41:49), today, discountsCap) - deltaC(i);
    sigma(47, i) = fzero(f, 0.2);
    tmp = vec(sigma(47, i));
    for k = 1:7
        sigma(47-k, i) = tmp(end-k);
    end
end

capPrev = capNext;

% Calibration from 12 to 15 years
capNext = priceCap(euribors(1:59), resetDates(1:61), strikes, vol(12, :), dates, discounts);
deltaC = capNext - capPrev;
discountsCap = interpolation(discounts, dates, resetDates(49:61));

for i = 1:m
    vec = @(sigmaTH) sigmaVec(sigmaTH, sigma, resetDates, 12, i, 12);
    % function to minimize
    f = @(sigmaTH) sumCaplet(euribors(48:59), strikes(i), vec(sigmaTH), resetDates(49:61), today, discountsCap) - deltaC(i);
    sigma(59, i) = fzero(f, 0.2);
    tmp = vec(sigma(59, i));
    for k = 1:11
        sigma(59-k, i) = tmp(end-k);
    end
end

capPrev = capNext;

% Calibration from 15 to 20 years
capNext = priceCap(euribors(1:end), resetDates(1:end), strikes, vol(13, :), dates, discounts);
deltaC = capNext - capPrev;
discountsCap = interpolation(discounts, dates, resetDates(61:end));

for i = 1:m
    vec = @(sigmaTH) sigmaVec(sigmaTH, sigma, resetDates, 13, i, 20);
    % function to minimize
    f = @(sigmaTH) sumCaplet(euribors(60:end), strikes(i), vec(sigmaTH), resetDates(61:end), today, discountsCap) - deltaC(i);
    sigma(end, i) = fzero(f, 0.2);

    tmp = interp1([resetDates(61),resetDates(end-1)],[sigma(59,i),sigma(end,i)],resetDates(61:end-1));
    for k = 1:19
        sigma(end-k, i) = tmp(end-k);
    end
end

end