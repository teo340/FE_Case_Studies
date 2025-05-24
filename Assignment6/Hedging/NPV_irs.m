function sensSwap = NPV_irs(dates, discounts, paymentDates, swapRate)
% NPV_irs: compute the DV01 of a payer swap to a 1 bp parallel rate shift
%
% INPUTS:
%   dates         - vector of dates defining the discount curve
%   discounts     - corresponding discount factors
%   paymentDates  - vector of fixed‑leg payment dates for the swap
%   swapRate      - fixed rate of the swap (par rate)
%
% OUPUT:
%   sensSwap      - sensitivity of the swap’s NPV to a 1 bp shift (DV01)

    % Define day count convention as European 30/360
    EU_30_360=6;

    % Compute year fractions for each fixed‑leg interval
    yf = yearfrac(paymentDates(1:end-1),paymentDates(2:end), EU_30_360);

    % Interpolate discount factors at each fixed‑leg payment date
    discounts_inter = interpolation(discounts,dates, paymentDates(2:end));

    % Calculate the Basis Point Value (BPV)
    BPV = sum(yf.*discounts_inter);

    % DV01 formula: change in NPV per 1 bp 
    % NPV(floating leg) ≈ 1 - discount factor at final payment
    % At inception the swap’s NPV is zero by definition (par swap)
    sensSwap = 1 - discounts_inter(end) - swapRate*BPV;
end