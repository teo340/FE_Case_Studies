function vec = sigmaVec(sigmaTH, sigma, resetDates, j, i, K)
% sigmaVec - Generates a vector of interpolated volatilities for a specific tenor
%
% This function creates a vector of volatilities by linear time interpolation
% between a previous known volatility and a target threshold volatility.
% It is used in calibration to generate consistent volatility structures
% for pricing caplets within a specific tenor period.
%
% Inputs:
%   sigmaTH - Target threshold volatility (for the end point of the tenor)
%   sigma - Current volatility matrix (partially filled during calibration)
%   resetDates - Reset dates for caplets
%   j - Current tenor index being processed
%   i - Current strike index being processed
%   K - Number of points to generate in the volatility vector
%
% Output:
%   vec - Vector of interpolated volatilities of length K

% Define day count convention (ACT_365)
ACT_365 = 3;

% Create interpolation function based on time fractions
% This performs linear interpolation in time between:
% - The last volatility from previous tenor: sigma(4*(j-1)-1,i)
% - The target threshold volatility: sigmaTH
sigma_inter = @(sigmaTH, target) sigma(4*(j-1)-1,i) + ...
              (sigmaTH - sigma(4*(j-1)-1,i)) / ...
              yearfrac(resetDates(4*(j-1)+1), resetDates(4*j+1), ACT_365) * ...
              yearfrac(resetDates(4*(j-1)+1), target, ACT_365);

% Initialize output vector
vec = zeros(K, 1);

% Fill vector with interpolated volatilities for intermediate points
for k = 1:K-1
    vec(k) = sigma_inter(sigmaTH, resetDates(4*(j-1)+k+1));
end

% Set last element to the target threshold volatility
vec(end) = sigmaTH;
end