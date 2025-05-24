function [x2, x5, x10] = hedging(sensIRS2, sensIRS5, sensIRS10, delta2y, delta5y, delta10y)
% hedgeCertificate computes the hedge ratios for mid market swaps  
% of different maturities (2y, 5y, 10y) to hedge the interest rate risk 
% (DV01 exposure) of a certificate.
%
% INPUTS:
%   sensIRS2   - Sensitivity of the 2-year IRS
%   sensIRS6   - Sensitivity of the 5-year IRS
%   sensIRS10  - Sensitivity of the 10-year IRS
%   delta2y    - DV01 contribution from the 2-year bucket
%   delta6y    - DV01 contribution from the 5-year bucket
%   delta10y   - DV01 contribution from the 10-year bucket
%   notional   - Notional amount of the certificate to be hedged
%
% OUTPUTS:
%   x2, x6, x10 - Notional amounts to trade in 2y, 5y, and 10y IRS respectively 
%                in order to hedge the DV01 exposure

notional = 100e6;
% Compute the hedge position in the 10-year IRS
% We offset the entire DV01 exposure using the 10y IRS first
x10 = -delta10y*notional / sensIRS10;

% Compute the hedge position in the 5-year IRS
% Offset the remaining exposure using the 5y IRS
x5 = -(sensIRS10 * x10 + delta5y*notional) / sensIRS5;


% Compute the hedge position in the 2-year IRS
% Offset the final remaining exposure using the 2y IRS
x2 = -(sensIRS10 * x10 + sensIRS5 * x5 + delta2y*notional) / sensIRS2;

end
