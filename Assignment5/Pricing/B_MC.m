function B_matrix = B_MC(dates, discounts, queryDates, vBMM, lambda, euribors, strikes)
% B_MC - Generates a Monte Carlo simulation of discount bond prices using the Bond Market Model
%
% This function performs Monte Carlo simulation to generate paths of discount bond prices
% following the Bond Market Model (BMM) dynamics. It incorporates correlation between
% consecutive time periods through the parameter lambda.
%
% Inputs:
%   dates - Dates for which discount factors are known
%   discounts - Known discount factors
%   queryDates - Dates for which to generate bond prices
%   vBMM - Calibrated BMM v-functions (by reset date and strike)
%   lambda - Mean reversion parameter for correlation
%   euribors - Forward Euribor rates
%   strikes - Strike rates used for volatility interpolation
%
% Output:
%   B_matrix - Matrix of simulated discount bond prices, dimensions (M x N)
%              where M is the number of simulations and N is the number of time points

% Set random seed for reproducibility
rng(42);

% Define day count convention (ACT_365)
ACT_365 = 3;

% Store valuation date
today = dates(1);

% Calculate time fractions from today to query dates
deltas = yearfrac(today, queryDates(3:end), ACT_365);

% Set number of Monte Carlo simulations
M = 1e7;

% Calculate year fractions between consecutive payment dates
yf = yearfrac(queryDates(2:end-2), queryDates(3:end-1), ACT_365);

% Determine number of time points for simulation
N = length(queryDates) - 2;

% Initialize matrix to store simulated bond prices
B_matrix = zeros(M, N);


% Loop through time periods to generate correlated bond price paths
for i = 1:N-1
    % Calculate correlation between consecutive periods using exponential decay
    rho = exp(-lambda * yf(i));
    
    % Generate correlated random normal variables
    Z = mvnrnd([0, 0], [1, rho; rho, 1], M);
    Z1 = Z(:, 1);
    Z2 = Z(:, 2);
    
    % Interpolate discount factors for relevant dates
    disc1 = interpolation(discounts, dates, queryDates(i+1));
    disc2 = interpolation(discounts, dates, queryDates(i+2));
    disc3 = interpolation(discounts, dates, queryDates(i+3));
    
    % Get strike rates and interpolate appropriate volatilities
    strike1 = euribors(i);
    v1 = spline(strikes, vBMM(i, :), strike1);
    
    strike2 = euribors(i+1);
    v2 = spline(strikes, vBMM(i+1, :), strike2);
    
    % Calculate bond price ratios using lognormal dynamics with correlation
    % B(t,T2)/B(t,T1) follows a geometric Brownian motion
    B_matrix(:, i) = disc2/disc1 * exp(-0.5*v1^2*deltas(i) + v1*sqrt(deltas(i))*Z1);
    B_matrix(:, i+1) = disc3/disc2 * exp(-0.5*v2^2*deltas(i+1) + v2*sqrt(deltas(i+1))*Z2);
end

% Add the initial discount bond B(0,3m) to complete the matrix
B_matrix = [interpolation(discounts, dates, queryDates(2))*ones(M, 1), B_matrix];
end