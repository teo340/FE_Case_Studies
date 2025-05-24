% Assignment 4 - GROUP 4
% AA 2024-2025
% runAssignment4_Group4

% clear workspace
clear all;
close all;
clc;

%% Settings and Data Import

formatData ='dd/MM/yyyy'; % Pay attention to your computer settings 

rng(42);   % Fix the random number generator ("the answer to Life, the Universe, and Everything")

load("ratesCurve.mat");
addpath("Bootstrap/");

% Bootstrap the input rate data to derive discount factors
[dates,discounts] = bootstrap(datesSet, ratesSet);


%% Point 1 - Certificate Spread Over Libor

% Define market parameters
rho = 0.45;
sigma1 = 0.162; sigma2 = 0.2;
d1 = 0.025; d2 = 0.029;
P = 0.9;
alpha = 1.1;
upfront = 0.03;
S0 = [100, 200];

% Set simulation start date and end date
startDate = dates(1); endDate = businessdayoffset(datetime("02-Feb-2027"));

% Set number of Monte Carlo simulations
M = 1e7;

% Generate simulated stock price paths using a Monte Carlo simulation function
[S1,S2] = stockPricesMC(rho, sigma1, sigma2, d1, d2, S0, startDate, discounts, dates, M);

% Define payment dates for Party B
partyB_dates = businessdayoffset(startDate+calyears(0:4));

% Interpolate discount factors at Party B payment dates
partyB_discounts = [1 interpolation(discounts, dates, partyB_dates(2:end))];

% Compute ratios (returns) between successive simulated stock prices for each asset
ratio1 = S1(:,2:end)./S1(:,1:end-1);
ratio2 = S2(:,2:end)./S2(:,1:end-1);

% Compute performance measure (St)
St = 0.25*sum(0.5*ratio1+0.5*ratio2, 2);

% Calculate the discounted payoff which is the excess performance over the protection level P (only positive differences)
disc_payoff = max(St-P,0);

% Compute the standard deviation of the discounted payoffs
stdev = std(disc_payoff);

% Set a confidence level
confLevel = 0.99;
% Compute confidence interval
IC = [mean(disc_payoff) - norminv(1-confLevel/2)*stdev,mean(disc_payoff) + norminv(1-confLevel/2)*stdev];

% Define payment dates for Party A
partyA_dates = businessdayoffset(startDate:calmonths(3):endDate);
% Interpolate discount factors for Party A payment dates
partyA_discounts = [1,interpolation(discounts,dates,partyA_dates(2:end))];
ACT_360 = 2;
partyA_deltas = yearfrac(partyA_dates(1:end-1),partyA_dates(2:end),ACT_360);
% Compute the Basis Point Value (BPV) by summing the products of discount factors and year fractions
BPV = sum(partyA_discounts(2:end).*partyA_deltas);

partyB_payment = alpha*mean(disc_payoff);
% Compute the spread over Libor (spol)
spol = (upfront+(partyB_payment+P)*partyB_discounts(end)-1)/BPV;
% Compute the confidence interval for the spread over Libor based on the payoff confidence interval
spol_IC = [(upfront+(alpha*IC(1)+P)*partyB_discounts(end)-1)/BPV, (upfront+(alpha*IC(2)+P)*partyB_discounts(end)-1)/BPV];

fprintf("Spread Over Libor: %d\n", spol);
fprintf("Confidence Interval(0.99): [%d, %d]\n",spol_IC(1),spol_IC(2));
fprintf("---------------\n");

%% Point 2 - Price of the Digital Option

notional = 1e7;

% spot and strike from data
load("cSelect.mat")
S_0 = cSelect.reference;
K = S_0;    % ATM
d = cSelect.dividend;

% Define the pricing date (t0) and maturity date, adjusting maturity for business days
t0 = datetime("02-Feb-2023");
maturity = datetime("02-Feb-2024");
maturity = businessdayoffset(maturity);

% Define the day-count convention
ACT_365 = 3;
T = yearfrac(t0, t0 + calyears(1), ACT_365);

% Interpolate the discount factor for one year from t0
discount_1y = interpolation(discounts, dates, t0 + calyears(1));

% Compute the forward price F0
F0  = S_0 / discount_1y * exp(-d * T);

% Obtain discount factor for option maturity using interpolation
discount_1y = interpolation(discounts, dates, maturity);

% Price the digital option using two methods:
% - 'black' represents the price computed using Black's formula,
% - 'smile' accounts for the volatility smile effect.
[black, smile] = digitalPrice(F0, K, maturity, notional, discount_1y);

fprintf("Digital Black: %d\n", black);
fprintf("Digital Price smile: %d\n", smile);
fprintf("---------------\n");

%% Point 3 - Mean Variance Mixture
% Define parameters for the mean-variance mixture model
sigma = 0.2;
eta = 3;
k = 1;
t = 1;

d = cSelect.dividend;

% Define a grid for log-moneyness values ranging from -25% to +25% (in steps of 1%)
logMoneyness = (-25:1:25) ./ 100;

% Compute call option prices using three different numerical integration methods
priceQuad = callPricing(logMoneyness, sigma, k, eta, 0, F0, discount_1y, t, 'Quad');
priceFFT = callPricing(logMoneyness, sigma, k, eta, 0, F0, discount_1y, t, 'FFT');
priceMC = callPricing(logMoneyness, sigma, k, eta, 0, F0, discount_1y, t, 'MC');

% Print the option price at ATM (where log moneyness equals zero)
mask = (logMoneyness == 0);
fprintf("------- Alpha = 0 -------\n");
fprintf("FFT: %d\n", priceFFT(mask));
fprintf("QUAD: %d\n",priceQuad(mask));
fprintf("MC: %d\n",priceMC(mask));
fprintf("-------------------------\n");

figure
plot(logMoneyness, priceQuad);
hold on; grid on;
plot(logMoneyness, priceFFT);
plot(logMoneyness, priceMC);
legend('Quad','FFT','MC');
hold off

%% Alpha = 1/3
alpha = 1/3;

% Recompute call prices using Quadrature and FFT methods with the new alpha value
priceQuad1 = callPricing(logMoneyness, sigma, k, eta, alpha, F0, discount_1y, t, 'Quad');
priceFFT1 = callPricing(logMoneyness, sigma, k, eta, alpha, F0, discount_1y, t, 'FFT');

% Print ATM prices for the alpha = 1/3 case for both methods
mask = (logMoneyness == 0);
fprintf("------- Alpha = 1/3 -------\n");
fprintf("FFT: %d\n", priceFFT1(mask));
fprintf("QUAD: %d\n",priceQuad1(mask));
fprintf("-------------------------\n");

figure
plot(logMoneyness, priceQuad1);
hold on; grid on;
plot(logMoneyness, priceFFT1);
legend('Quad','FFT');
hold off

% Compute the average difference between the FFT prices of the two model variants
dist = mean(priceFFT1 - priceFFT); 
relDiff = dist/priceFFT(mask);
disp(relDiff)

%% Point 4 - Volatility Surface Calibration

strikes = cSelect.strikes;
vol = cSelect.surface;

% Calculate market option prices using Black's formula for each strike,
% discounting by the interpolated factor
realPrices = discount_1y*blkprice(F0, strikes, 0, t, vol);

% Compute log-moneyness for each strike (log of the ratio of F0 to strike)
logMoneyness = log(F0./strikes);

% Set the model parameter alpha to 2/3 for the calibration step
alpha = 2/3;

% Define a function handle 'prices' which returns model call prices
% computed via FFT given a vector p = [sigma, k, eta] of parameters
prices = @(p) callPricing(logMoneyness, p(1), p(2), p(3), alpha, F0, discount_1y, t, 'FFT');

%% Point 4 - Distance minimization

% create the distance function to minimize
dist = @(x) sum((prices(x) - realPrices).^2);

% create the constraint
const = @(x) constraint(x, alpha);

% calibrate the model using fmincon
% initial guess
x0 = [0.2, 1, 3];

% Set up linear constraints to ensure parameters remain positive
A = [
    -1, 0, 0;
    0, -1, 0;
    0, 0, -1;
];
b = [
    0;
    0;
    0;
];

% Calibrate the model parameters using constrained nonlinear optimization (fmincon)
options = optimoptions('fmincon', 'Display', 'off');
[x, fval] = fmincon(dist, x0, A, b, [], [], [], [], [], options);

%% Point 4 - Implied Volatilities

% Compute model-generated option prices using the calibrated parameters
modelPrices = prices(x);
% Calculate the implied volatility using the Black model
impvol = blkimpv(F0, strikes, 0, t, modelPrices./discount_1y);

plot(strikes, vol);
hold on; grid on;
plot(strikes, impvol);
legend("Black", "Exp-Levy");
title("Implied volatility comparison");