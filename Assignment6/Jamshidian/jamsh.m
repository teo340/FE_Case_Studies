function [price] = jamsh(dates, discounts, dates_years, coupon, sigma, a)
    
ACT_360 = 2;
ACT_365 = 3;

today = dates(1);
t_alpha = dates_years(1);
t_omega = dates_years(end);

deltas = yearfrac(dates_years(1:end-1), dates_years(2:end), ACT_360);

yf = yearfrac(today, dates_years, ACT_365);

coupons = coupon * deltas;
coupons(end) = 1 + coupons(end);
disc = interpolation(discounts, dates, dates_years);
fwd_disc = disc(2:end)/disc(1);
sigmaHJM = @(s,T) sigma/a * (1 - exp(-a * (T - s)));

fwd_stoc_disc = cell(1,length(coupons));

for i = 1:length(coupons)
    sigma_prev = sigmaHJM(yf(1),yf(1));
    sigma_next = sigmaHJM(yf(1),yf(i+1));
    integrand =  @(t) sigmaHJM(t,yf(1)).^2 - sigmaHJM(t,yf(i+1)).^2;
    intSigma = quadgk(integrand, 0, yf(1));
    fwd_stoc_disc{i} = @(x) fwd_disc(i) * exp(-x/sigma * (sigma_next-sigma_prev)-0.5*intSigma);
end


price_CB = @(x) - 1;
for i = 1:length(coupons)
    price_CB = @(x) price_CB(x) + fwd_stoc_disc{i}(x)*coupons(i);
end

x_star = fzero(price_CB,0);

K_i = zeros(1,length(coupons));
for i = 1:length(coupons)
    sigma_prev = sigmaHJM(yf(1),yf(1));
    sigma_next = sigmaHJM(yf(1),yf(i+1));
    integrand =  @(t) sigmaHJM(t,yf(1)).^2 - sigmaHJM(t,yf(i+1)).^2;
    intSigma = quadgk(integrand, 0, yf(1));
    K_i(i) = fwd_disc(i) * exp(-x_star/sigma * (sigma_next-sigma_prev) - 0.5*intSigma);
end

Calls = zeros(1,length(K_i));
Puts = Calls;
for i = 1:length(K_i)
    Calls(i) = priceCall(K_i(i), fwd_disc(i), disc(1), today, t_alpha, dates_years(i+1), sigma, a);
    Puts(i) = Calls(i) - disc(1)*(fwd_disc(i)-K_i(i));
end

price = sum(coupons.*Puts);

end
