% runAssignment2
%  group X, AY20ZZ-20ZZ
% Computes Euribor 3m bootstrap with a single-curve model
%
% This is just code structure: it should be completed & modified (TBM)
%
% to run:
% > runAssignment2_TBM

clear all;
close all;
clc;

%% Settings
formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

%% Read market data
% This fuction works on Windows OS. Pay attention on other OS.

%[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap.xls', formatData);

%% Bootstrap
% dates includes SettlementDate as first date
[dates, discounts]=bootstrap(datesSet, ratesSet);

% the df at 1 year it's not supposed to be in the bootstrap,
% yet it's needed later (for instance to bootstrap swaps).

%% Compute Zero Rates
zRates = 100*zeroRates(dates, discounts);

%% Plot Results

%discount curve
%dates = sort(dates);
figure;
plot(dates, discounts,'b-*', 'LineWidth', 1.5);
grid on;
title('Discount Factors');

%zero-rates
figure;
plot(dates,zRates,'b-*', 'LineWidth', 1.5); 
grid on;
title('Zero Rates');

%% Pricing IB Bond
% the parameter df_1y has been saved above, however
% it's not included in the bootstrapped curve.
bondprice = price_bond(datesSet, ratesSet);
fprintf('Price of the bond: %0.15f\n', bondprice);
fprintf('----------\n');

%% Sensititivities
today = datetime('02-Feb-2023', 'InputFormat', 'dd-MMM-yyyy');

rates_shifted = rates_shift(ratesSet);
[~, df_shifted] = bootstrap(datesSet, rates_shifted);

fixedLeg = [datetime('2-Feb-2024'); dates(13:18)];
fixedRate = 0.028175;
[DV01,BPV,DV01_z] = sensSwap(today, fixedLeg, fixedRate, dates', discounts, df_shifted);

MacD = sensCouponBond(today, fixedLeg, fixedRate, dates, discounts);

fprintf('DV01: %d\n',DV01);
fprintf('BPV: %d\n',BPV);
fprintf('DV01_z: %d\n',DV01_z);
fprintf('MacD: %d\n', MacD);
fprintf('----------\n');

%% Exercise 4
NPV_1 = compute_NPV(discounts, dates, 1500, 0.05);
NPV_2 = compute_NPV(discounts, dates, 6000, 0.05);

fprintf('Initial CF: 1500 -> %0.6f\n', NPV_1);
fprintf('Initial CF: 6000 -> %0.6f\n', NPV_2);
fprintf('----------\n');