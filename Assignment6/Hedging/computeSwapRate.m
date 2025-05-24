
function swapRate = computeSwapRate(dates, discounts, paymentDates)
% ComputeSwapRate Calculates the swap rate for an interest rate swap
%   swapRate = computeSwapRate(dates, discounts, paymentDates) computes the 
%   fixed rate that makes the net present value of the swap equal to zero
%
%   INPUTS:
%       dates - Vector of dates for discount factors (date format)
%       discounts - Vector of discount factors corresponding to dates
%       paymentDates - Vector of dates when swap payments occur (date format)
%
%   OUTPUT:
%       swapRate - The calculated fixed swap rate

    % Define day count convention as European 30/360
    EU_30_360=6;

    % Calculate year fractions between consecutive payment dates
    yf = yearfrac(paymentDates(1:end-1),paymentDates(2:end), EU_30_360);

    % Interpolate discount factors for payment dates
    discounts_inter = interpolation(discounts,dates, paymentDates(2:end));

    % Calculate Basis Point Value (BPV)
    BPV = sum(yf.*discounts_inter);

    % Calculate swap rate as (1 - final discount factor) / BPV   
    swapRate = (1-discounts_inter(end))/BPV;
end