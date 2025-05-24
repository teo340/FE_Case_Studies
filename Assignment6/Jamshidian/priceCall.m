function price = priceCall(strike, fwd, disc, t0, start, expiry, sigma, a)

ACT_365 = 3;

% compute the volatility using the integral
talpha = yearfrac(t0, start, ACT_365);
ttm = yearfrac(t0, expiry, ACT_365);

% compute the volatility integral
sigmaHJM = @(s,T) sigma/a * (1 - exp(-a * (T - s)));
v_2 = @(t) (sigmaHJM(t, ttm) - sigmaHJM(t,talpha)).^2;

% compute the integral
V = sqrt(1 / talpha * quadgk(v_2, 0, talpha));

% define d1 and d2 for every T_i given
d1 = log(fwd / strike) / (V * sqrt(talpha)) + 0.5 * V * sqrt(talpha);
d2 = d1 - V * sqrt(talpha);
price = disc * (fwd * normcdf(d1) - strike * normcdf(d2));

end