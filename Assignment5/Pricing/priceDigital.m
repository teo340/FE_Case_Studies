function price = priceDigital(euribor, mkt_strikes, strike, vols, t0, tprev, tnext, B, flag)
% PriceDigital calculates the price of a digital option on interest rates
%   price = priceDigital(euribor, mkt_strikes, strike, vols, t0, tprev, tnext, B, flag)
%   computes the price of a digital option using either the standard Black model
%   or with smile adjustment
%
%   INPUTS:
%       euribor - Forward Euribor rate
%       mkt_strikes - Vector of market strike for volatility surface
%       strike - Specific strike rate for the digital option
%       vols - Volatilities corresponding to mkt_strikes
%       t0 - Valuation date
%       tprev - Reset date of the option
%       tnext - Payment date of the option
%       B - Discount factor for the payment date
%       flag - Pricing method: "Black" for standard model or "DigitalRisk" for risk-adjusted
%
%   OUTPUT:
%       price - Price of the digital option

% Day‑count convention for volatility (Act/365)
ACT_365 = 3;

% Calculate time to reset date using Actual/365
deltaT = yearfrac(t0, tprev, ACT_365);

% Interpolate volatility at the exact strike level
int_vol = spline(mkt_strikes, vols, strike);

% Black d2 term
d2 = log(euribor/strike)/(int_vol*sqrt(deltaT))-0.5*int_vol*sqrt(deltaT);


% Standard Black digital price = discount × N(d2)
if flag == "Black"
    price = B*normcdf(d2);
% Risk-adjusted pricing that accounts for volatility skew
elseif flag == "DigitalRisk"
    % Locate segment on vol surface containing the strike
    idx = find(mkt_strikes(1:end-1) <= strike  & mkt_strikes(2:end) >= strike, 1);
    if isempty(idx)
        error('Strike not found in the vol surface');
    end

    % Estimate local slope of the smile: Δσ/ΔK between adjacent strikes
    sigma_diff = (vols(idx+1)-vols(idx))/(mkt_strikes(idx+1)-mkt_strikes(idx));

    % Adjusted digital price = Black price minus smile‑risk term
    % Smile‑risk term = slope × caplet vega
    price = B*normcdf(d2) - sigma_diff*vegaCaplet(euribor, strike, int_vol, t0, tprev, tnext, B);
else
    error("Flag not defined");
end

end