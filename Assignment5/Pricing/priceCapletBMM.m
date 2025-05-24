function price = priceCapletBMM(euribor, strikes, v, t0, tprev, tnext, B)
% priceCapletBMM: price a caplet under the Bond Market Model (BMM)
%
% INPUTS:
%       euribor - forward Euribor rate
%       strikes - vector of strike rates for the caplet
%       v - BMM volatility parameter 
%       t0 - valuation date
%       tprev - reset date of the caplet
%       tnext - payment date of the caplet
%       B - discount factor for the payment date
%
% OUTPUT:
%   price     - vector of caplet prices for each strike under BMM

% Define day count conventions
ACT_360 = 2;  % Actual/360 day count convention for interest rate calculation
ACT_365 = 3;  % Actual/365 day count convention for year fraction in volatility scaling

% Calculate year fraction between reset and payment dates using Actual/360
yf = yearfrac(tprev, tnext, ACT_360);

% Calculate time to reset date in years using Actual/365
deltaT = yearfrac(t0, tprev, ACT_365);

% Calculate d1 parameter for BMM
d1= log((1+yf*euribor)./(1+yf*strikes))./(v*sqrt(deltaT)) + 0.5*v*sqrt(deltaT);
% Calculate d2 parameter for BMM
d2= log((1+yf*euribor)./(1+yf*strikes))./(v*sqrt(deltaT)) - 0.5*v*sqrt(deltaT);

% Caplet price under BMM
price = B*((1+yf*euribor).*normcdf(d1)-(1+yf.*strikes).*normcdf(d2));

end