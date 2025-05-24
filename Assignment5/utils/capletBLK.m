function caplet_prices = capletBLK(mkt_strike, mkt_vols, today, dates, discounts)
% This function compute the prices of a list of caplet, using Black Formula.
% It's used for the sole porpose of calibration, NOT for pricing, for which instead
% the function Pricing/priceCaplet.m is used.#
% INPUTS:
% mkt_strike: strikes at which to price the caplets
% mkt_vols: volatility to price the caplets
% today: settlement date, datetime
% dates: bootstrap dates, datetime
% discounts: bootstrap discounts

% compute the Euribor
discounts_fwd = discounts(2:end) ./ discounts(1:end-1);
year_frac_from_t0 = yearfrac(today, dates(1:end-1), 3);
year_frac_moving = yearfrac(dates(1:end-1), dates(2:end), 2);
L_fwd = (1 ./ discounts_fwd - 1) ./ year_frac_moving;

% change dimension for blkprice compatibility
L_fwd = L_fwd';
year_frac_from_t0 = year_frac_from_t0';
year_frac_moving = year_frac_moving';
discounts = discounts';

% use blkprice to quote the instrument
caplet_prices = transpose(blkprice(L_fwd, mkt_strike, 0, year_frac_from_t0,mkt_vols).* year_frac_moving .* discounts(2:end));

end
