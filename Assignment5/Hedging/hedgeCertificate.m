function [x2, x6, x10] = hedgeCertificate(sensIRS2, sensIRS6, sensIRS10, delta2y, delta6y, delta10y, notional)
% hedgeCertificate computes the hedge ratios for interest rate swaps (IRS) 
% of different maturities (2y, 6y, 10y) to hedge the interest rate risk 
% (DV01 exposure) of a certificate.
%
% INPUTS:
%   sensIRS2   - Sensitivity of the 2-year IRS
%   sensIRS6   - Sensitivity of the 6-year IRS
%   sensIRS10  - Sensitivity of the 10-year IRS
%   delta2y    - DV01 contribution from the 2-year bucket
%   delta6y    - DV01 contribution from the 6-year bucket
%   delta10y   - DV01 contribution from the 10-year bucket
%   notional   - Notional amount of the certificate to be hedged
%
% OUTPUTS:
%   x2, x6, x10 - Notional amounts to trade in 2y, 6y, and 10y IRS respectively 
%                in order to hedge the DV01 exposure

% Total DV01 exposure for the 10-year maturity bucket
DV01_10y = delta2y + delta6y + delta10y;

% Compute the hedge position in the 10-year IRS
% We offset the entire DV01 exposure using the 10y IRS first
x10 = -DV01_10y * notional / sensIRS10;

% Residual DV01 (after hedging with x10) for the 6-year bucket
DV01_6y = delta2y + delta6y;

% Compute the hedge position in the 6-year IRS
% Offset the remaining exposure using the 6y IRS
x6 = -(sensIRS10 * x10 + DV01_6y * notional) / sensIRS6;

% Residual DV01 for the 2-year bucket
DV01_2y = delta2y;

% Compute the hedge position in the 2-year IRS
% Offset the final remaining exposure using the 2y IRS
x2 = -(sensIRS10 * x10 + sensIRS6 * x6 + DV01_2y * notional) / sensIRS2;

end
