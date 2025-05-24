function [blackPrice, smilePrice] = digitalPrice(F0, strike, maturity, notional, B)
% This function computes the price of a digital call option, using
% two different models: black and smile-adjusted.
% INPUTS:
% F0: forward price
% strike: strike price
% maturity: maturity date
% notional: notional amount
% B: discount factor
% OUTPUTS:
% blackPrice: price using the Black model
% smilePrice: price using the smile-adjusted model

% Vol smile information is loaded from a .mat file
load("cSelect.mat")

strikes = cSelect.strikes;
vol = cSelect.surface;
today = datetime("02-Feb-2023");

payment = 0.07*notional;

ACT_365 = 3;
deltaT = yearfrac(today, maturity, ACT_365);

% implied volatility as a spline interpolation of the smile
impvol = @(K) interp1(strikes, vol, K, 'spline');

d1 = @(K) (log(F0 / K) + (0.5 * impvol(K)^2) * deltaT) / (impvol(K) * sqrt(deltaT));
d2 = @(K) d1(K) - impvol(K) * sqrt(deltaT);

blackPrice = B*payment*normcdf(d2(strike));
% fprintf("Black Price: %d\n", blackPrice);

% Vega of the call option
vega = @(K) F0*exp(-d1(K).^2/2)/sqrt(2*pi)*sqrt(deltaT);

% We take the derivative of the smile, approximating it with a finite difference.
% We need to find the two closest strikes to the one we are interested in.
i = find(strikes(1:end-1) <= strike  & strikes(2:end) >= strike, 1);
if isempty(i)
    error('Strike not found in the vol surface');
end
sigma_diff = (vol(i+1)-vol(i))/(strikes(i+1)-strikes(i));

smilePrice = B*payment*(normcdf(d2(strike))-sigma_diff*vega(strike));
% fprintf("Smile Price: %d\n", smilePrice*B*payment);

end