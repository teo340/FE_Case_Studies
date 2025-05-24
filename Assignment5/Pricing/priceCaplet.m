function price = priceCaplet(euribor, strike, vol, t0, tprev, tnext, B)
% priceCaplet: price a single caplet under the Libor Market Model (LMM)
%
% INPUTS:
%   euribor   - forward Euribor rate
%   strike    - vector of strike rates for the caplet
%   vol       - vector of volatilities corresponding to each strike 
%   t0        - valuation date
%   tprev     - reset date of the caplet
%   tnext     - payment date of the caplet
%   B         - discount factor to the payment date
%
% OUTPUT:
%   price     - vector of caplet prices for each strike


% Define day count conventions
ACT_360 = 2;  % Actual/360 day count convention for interest rate calculation
ACT_365 = 3;  % Actual/365 day count convention for volatility timing

% Calculate year fraction between reset and payment dates using Actual/360
yf = yearfrac(tprev, tnext, ACT_360);

% Calculate time to reset date in years using Actual/365
deltaT = yearfrac(t0, tprev, ACT_365);

% Calculate d1 parameter for LMM
d1 = log(euribor./strike)./(vol*sqrt(deltaT))+0.5*vol*sqrt(deltaT);
% Calculate d2 parameter for LMM
d2 = log(euribor./strike)./(vol*sqrt(deltaT))-0.5*vol*sqrt(deltaT);

% LMM caplet price
price = yf*B*(euribor.*normcdf(d1)-strike.*normcdf(d2));

end