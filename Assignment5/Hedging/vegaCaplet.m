function vega = vegaCaplet(euribor, strike, vol, t0, tprev, tnext, B)
% vegaCaplet computes the Black vega of a caplet, i.e., the sensitivity of 
% the caplet price to changes in volatility.
%
% INPUTS:
%   euribor  - Value of Euribor 
%   strike   - Strike rate of the caplet
%   vol      - Volatility 
%   t0       - Valuation date
%   tprev    - Start date of the caplet
%   tnext    - End date of the caplet
%   B        - Discount factor 
%
% OUTPUT:
%   vega     - Vega of the caplet (sensitivity of price w.r.t. volatility)

% Define day count conventions
ACT_360 = 2;  
ACT_365 = 3;  

% Compute the year fraction for the caplet period using ACT/360
yf = yearfrac(tprev, tnext, ACT_360);  

% Compute the time to expiry from the valuation date to the caplet start date
deltaT = yearfrac(t0, tprev, ACT_365);  
% Compute the d1 term from the Black formula
d1 = log(euribor / strike) / (vol * sqrt(deltaT)) + 0.5 * vol * sqrt(deltaT);

% Calculate vega using Black's formula for caplets
vega = euribor * yf * B * sqrt(deltaT) * exp(-d1^2 / 2) / sqrt(2 * pi);

end
