function [x6, x10] = hedgeVega(vega6y, vega10y, vega_cap6y, vega_cap10y, notional)
% hedgeVega computes the hedge notional amounts needed in 6-year and 10-year  ATM-Swap Rate Caps
% to neutralize the vega exposure.
%
% INPUTS:
%   vega6y        - Vega exposure in the 6-year bucket
%   vega10y       - Vega exposure in the 10-year bucket
%   vega_cap6y    - Vega sensitivity of a 6-year cap 
%   vega_cap10y   - Vega sensitivity of a 10-year cap 
%   notional      - Notional amount to be hedged
%
% OUTPUTS:
%   x6, x10       - Notional amounts to trade in 6y and 10y ATM-Swap Rate caps, respectively, 
%                   in order to hedge the product's vega exposure

% Step 1: Hedge the 10-year vega exposure using the 10y cap
% Compute the notional of 10y caps needed to offset the 10y vega exposure
x10 = -notional * vega10y / vega_cap10y;

% Step 2: Hedge the residual 6-year vega exposure using the 6y cap
% The residual vega includes the original 6y exposure and the indirect impact
% of the x10 hedge on the 6y bucket 
x6 = -(x10 * vega_cap10y + vega6y * notional) / vega_cap6y;

end
