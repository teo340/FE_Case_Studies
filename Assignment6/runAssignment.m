clc;
clear;
close all;

formatData ='dd/MM/yyyy'; % Pay attention to your computer settings 
addpath("./Bootstrap");
addpath("./Pricing/");
addpath("./Swaption/");
addpath("./Hedging/");
addpath("./Jamshidian/");

%% Data and Bootstrap
load("Dataset/ratesCurve.mat");
load("Dataset/cSelect.mat");

[dates,discounts] = bootstrap(datesSet, ratesSet);

%% Computing Upfront - 2year certificate
spol = 0.015;
upfront_VG = computeUpfront(spol, dates, discounts, "VG");
upfront_BLK = computeUpfront(spol, dates, discounts, "Black");
upfront_Smile = computeUpfront(spol, dates, discounts, "Smile");

fprintf("Upfront VG: %d\n", upfront_VG);
fprintf("Upfront BLK: %d\n", upfront_BLK);
fprintf("Upfront Smile: %d\n", upfront_Smile);

%% Computing Upfront - 3year certificate
% N = 1e5
[upfront, CI] = computeUpfront3y(0.015, dates, discounts);

fprintf("Upfront 3y MC: %d\n", upfront);
fprintf("Confidence Interval: [%d, %d]\n", CI(1), CI(2));

%% Swaption Pricing - Hull White Tree Approach

a = 10/100;
sigma = 0.8/100;

strike = 0.05;
deltaT = 10/365;
t_alpha = 2; t_omega = 10;

tic;
[price, tree] = swaptionPrice(t_omega, strike, dates, discounts, deltaT);
toc;

% display the value of the swaption
fprintf("--- Bermudan swaption value ---\n")
fprintf("The value of the Bermudan swaption is: %d\n", price);

%% Swaption Pricing - Convergence Check

N = [40, 30, 20, 10, 5, 1];

% Warning! Using a daily time step 1/365 may take a while to run, around 10
% seconds. To avoid performance issues, we suggest to run the following
% script to to N = 5, where the convergence is already quite visible.

deltaT = N/365;
prices = zeros(1,length(deltaT));
for i=1:length(deltaT)
    [prices(i), tree] = swaptionPrice(t_omega, strike, dates, discounts, deltaT(i));
    fprintf("deltaT = %d/365 ---> Swaption = %d\n", N(i), prices(i));
end
%% Jamshidian for lower and upper bound of bermudan swaption
clc
enddate = datetime('02-Feb-2033');
dates_years = businessdayoffset(dates(1):calyears(1):enddate);
price = zeros(1, 8);
for i = 2:9
    price(i-1) = jamsh(dates, discounts, dates_years(i+1:end), strike, sigma, a);
end

lower_bound = max(price);
upper_bound = sum(price);

fprintf('lower_bound: %.8f\n', lower_bound)
fprintf('upper_bound: %.8f\n', upper_bound)

%% Buckets - Rates Shift
deltaT=10/365;
[price, tree] = swaptionPrice(t_omega, strike, dates, discounts, deltaT);

ratesShifted2y = bucket2y(ratesSet);
ratesShifted5y = bucket5y(ratesSet);
ratesShifted10y = bucket10y(ratesSet);

[~, discounts2y] = bootstrap(datesSet, ratesShifted2y);
[~, discounts5y] = bootstrap(datesSet, ratesShifted5y);
[~, discounts10y] = bootstrap(datesSet, ratesShifted10y);

%% Hedging - Sensitivities

[price2y, tree2y] = swaptionPrice(t_omega, strike, dates, discounts2y, deltaT);
[price5y, tree5y] = swaptionPrice(t_omega, strike, dates, discounts5y, deltaT);
[price10y, tree10y] = swaptionPrice(t_omega, strike, dates, discounts10y, deltaT);

deltaPrice2y = price2y - price;
deltaPrice5y = price5y - price;
deltaPrice10y = price10y - price;

%% Hedging - IRS

resetDates = dates(1):calyears(1):dates(1)+calyears(t_omega);
resetDates = businessdayoffset(resetDates);
swapRate2y= computeSwapRate(dates, discounts, resetDates(1:2));
swapRate5y= computeSwapRate(dates, discounts, resetDates(1:5));
swapRate10y= computeSwapRate(dates, discounts, resetDates(1:end));

ratesSetShifted = ratesShiftAll(ratesSet);
[~, discountsShifted] = bootstrap(datesSet, ratesSetShifted);
sensSwap2y = NPV_irs(dates, discountsShifted, resetDates(1:2) ,swapRate2y);
sensSwap5y = NPV_irs(dates, discountsShifted, resetDates(1:5) ,swapRate5y);
sensSwap10y = NPV_irs(dates, discountsShifted, resetDates(1:end),swapRate10y);

[x2, x5, x10] = hedging(sensSwap2y, sensSwap5y, sensSwap10y, deltaPrice2y, deltaPrice5y, deltaPrice10y);

fprintf("IRS 2y: %d\n", x2);
fprintf("IRS 5y: %d\n", x5);
fprintf("IRS 10y: %d\n", x10);
