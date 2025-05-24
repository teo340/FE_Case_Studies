function upfront = computeUpfront(euribor, strikes, volSurface, spol, resetDates, dates, discounts, flag)
% ComputeUpfront calculates the upfront for a structured bond
%
%   INPUTS:
%       euribor - Vector of Euribor rates
%       strikes - Vector of strike prices for volatility surface
%       volSurface - Matrix of volatilities indexed by time and strike
%       spol - Spread value for payment calculation
%       resetDates - Vector of reset dates for the instrument
%       dates - Vector of dates for discount factors
%       discounts - Vector of discount factors corresponding to dates
%       flag - Flag for digital option pricing method
%
%   OUTPUT:
%       upfront - The calculated upfront payment amount


% Store today's date as the first date in the dates vector
today = dates(1);

% Interpolate discount factors for reset dates
discounts_inter = interpolation(discounts, dates, resetDates(2:end));

% Define day count convention as Actual/360
ACT_360 = 2;
% Calculate year fractions between consecutive reset dates
yf = yearfrac(resetDates(1:end-1), resetDates(2:end), ACT_360);

% Calculate Basis Point Value (BPV)
BPV = sum(yf.*discounts_inter);

% Calculate payment for Bank XX
partyXXpayment = 1-discounts_inter(end)+spol*BPV;

% Define cap strike prices for different periods
Kcap1 = 0.043;  % Cap strike for first period
Kcap2 = 0.046;  % Cap strike for second period
Kcap3 = 0.052;  % Cap strike for third period

% Define digital strike prices for different periods
Kdig1 = 0.043;  % Digital strike for first period
Kdig2 = 0.046;  % Digital strike for second period
Kdig3 = 0.052;  % Digital strike for third period

% Get the number of reset dates
N = length(resetDates);

% Initialize coupon vector to store calculated coupon values
coupon = zeros(1,N-1);

% Set first coupon as a fixed 3% rate adjusted by discount factor
coupon(1) = discounts_inter(1)*0.03;

% Calculate coupons for first period (up to (and including) the 3rd year)
for i=2:12
    % Interpolate discount factor for payment date
    B = interpolation(discounts, dates, resetDates(i+2));

    % Price digital option
    digital = priceDigital(euribor(i), strikes, Kdig1, volSurface(i,:), today, resetDates(i+1), resetDates(i+2), B, flag);

    % Interpolate volatility for the cap strike Kcap1
    int_vol = spline(strikes, volSurface(i,:), Kcap1);

    % Price cap option
    cap = priceCap(euribor(1:i), resetDates(1:i+2), Kcap1, int_vol, dates, discounts);

    % Calculate year fraction for this period
    yf = yearfrac(resetDates(i+1), resetDates(i+2), ACT_360);

    % Calculate coupon
    coupon(i) = yf*discounts_inter(i)*0.011+(discounts_inter(i)-discounts_inter(i+1))-cap-0.009*digital;
end

% Calculate coupons for second period (after 3y and up to (and including) the 6y)
for i=13:24
   % Interpolate volatility for the cap strike Kcap2
    int_vol = spline(strikes, volSurface(i,:), Kdig2);
    
    % Interpolate discount factor for payment date
    B = interpolation(discounts, dates, resetDates(i+2));
    
    % Price digital option using strike Kdig2
    digital = priceDigital(euribor(i), strikes, Kdig2, volSurface(i,:), today, resetDates(i+1), resetDates(i+2), B, flag);
    
    % Price cap option
    cap = priceCap(euribor(1:i), resetDates(1:i+2), Kcap2, int_vol, dates, discounts);
    
    % Calculate year fraction for this period
    yf = yearfrac(resetDates(i+1), resetDates(i+2), ACT_360);
    
    % Calculate coupon
    coupon(i) = yf*discounts_inter(i)*0.011+(discounts_inter(i)-discounts_inter(i+1))-cap-0.009*digital;
end

% Add an additional reset date (May 2, 2033) for the final calculations
resetDates = [resetDates, businessdayoffset(datetime("02-May-2033"))];

% Interpolate discount factor for the additional reset date
discounts_inter = [discounts_inter, interpolation(discounts, dates, resetDates(end))];

% Calculate coupons for third period (after 6y and up to 10y)
for i=25:40
    % Interpolate volatility for the cap strike Kcap3
    int_vol = spline(strikes, volSurface(i,:), Kdig3);
    
    % Interpolate discount factor for payment date
    B = interpolation(discounts, dates, resetDates(i+2));
    
    % Price digital option using strike Kdig3
    digital = priceDigital(euribor(i), strikes, Kdig3, volSurface(i,:), today, resetDates(i+1), resetDates(i+2), B, flag);
    
    % Price cap option
    cap = priceCap(euribor(1:i), resetDates(1:i+2), Kcap3, int_vol, dates, discounts);
    
    % Calculate year fraction for this period
    yf = yearfrac(resetDates(i+1), resetDates(i+2), ACT_360);
    
    % Calculate coupon: fixed rate + forward rate - cap adjustment - digital adjustment
    coupon(i) = yf*discounts_inter(i)*0.011+(discounts_inter(i)-discounts_inter(i+1))-cap-0.009*digital;
end

% Calculate upfront payment as the difference between Bank XX payment and sum of all coupons
upfront = partyXXpayment - sum(coupon);
end