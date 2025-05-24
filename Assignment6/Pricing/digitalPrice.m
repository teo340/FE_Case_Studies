function [price] = digitalPrice(strike, sigma, k, eta, alpha, F0, B, T)

% We price a digital as an infinitesimal bull spread

epsilon = 1e-3;
call1 = callPricing(log(F0/strike), sigma, k, eta, alpha, F0, B, T, 'FFT');
call2 = callPricing(log(F0/(strike+epsilon)), sigma, k, eta, alpha, F0, B, T, 'FFT');

price = 1/epsilon*(call1-call2);

end