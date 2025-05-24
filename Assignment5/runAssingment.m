clc;
clear;
close all;

addpath("Bootstrap/");
addpath("Pricing/");
addpath("Calibration");
addpath("Hedging/");
addpath("utils/");

%% Loading the relevant data

load("Dataset/ratesCurve.mat");
load("Dataset/mktData.mat");

[dates, discounts] = bootstrap(datesSet, ratesSet);
mkt_strikes = mktData.strikes;
mkt_years = mktData.years;
mkt_vols = mktData.vols;

%% Point 1.1 - Calibrating the Libor Market Model (LMM)

% Our structured bond is based on the 3M Euribor. We set up the reset dates,
% starting from today, 02-Feb-2023, up to twenty years in the future.

today = dates(1);
resetDates = today:calmonths(3):datetime("02-Feb-2043");
resetDates = businessdayoffset(resetDates);

discountsReset = interpolation(discounts, dates, resetDates);
discountsReset(1) = 1;

ACT_360 = 2;
fwd_discounts = discountsReset(2:end)./discountsReset(1:end-1);
deltas = yearfrac(resetDates(1:end-1),resetDates(2:end),ACT_360);

euribor = 1./deltas.*(1./fwd_discounts-1);
euribor = euribor(2:end);

%% Point 1.1 - Calibrating the Libor Market Model (LMM)

% The calibration to quoted market prices of caps takes around 2.5 seconds.
% Once the fitting process is complete, we plot the obtained volatility
% surface and compare it with the flat mkt_vol. The two are quite similar,
% but the former is much more granular.

tic;
sigma = calibrateLMM(euribor, mkt_strikes, resetDates, mkt_vols, dates, discounts);
toc;

figure;
[T,K] = meshgrid(resetDates(2:end-1), mkt_strikes);
surf(T,K,sigma')
title("Libor Market Model - Volatility Surface");

figure;
[T1,K1] = meshgrid(mkt_years, mkt_strikes);
surf(T1,K1,mkt_vols');
title("Market Data - Flat Vol");

%% Point 1.2 - Determining the upfront

% To price the contract, we don't need the full data until 20 years, but
% only up to the first 10. So we reduce the vol matrix and the rates.

euribor = euribor(1:40);
resetDates = resetDates(1:41);
sigma = sigma(1:40,:);

% Given a spread over libor of 2%, we compute the upfront of the contract
% ignoring and consindering the Digital Risk

spol = 0.02;
upfront_black = computeUpfront(euribor, mkt_strikes, sigma, spol, resetDates, dates, discounts, "Black");
upfront_smile = computeUpfront(euribor, mkt_strikes, sigma, spol, resetDates, dates, discounts, "DigitalRisk");

fprintf("Upfront Black: %d\n", upfront_black);
fprintf("Upfront smile: %d\n", upfront_smile);

%% Point 1.3 - Delta Bucket Sensitivities

% For each instrument of our bootstrapped curve, we study the sensitivity
% of our contract. We manually bump each rate by 1bp and compute the
% difference in terms of upfront.

% Warning! Due to the high number of complex operations, like bootstrapping
% and volatility calibration, this block may take some minutes to run.

tic;
imax=find(dates == '02-Feb-2033');
upfrontDiff=zeros(1,imax);
resetDates_all = today:calmonths(3):datetime("02-Feb-2043");
resetDates_all = businessdayoffset(resetDates_all);

for i=1:imax
    fprintf("%0.2f done\n", i/imax*100);

    shiftedRates = ratesShift(ratesSet, datesSet, datenum(dates(i+1)), 1e-4);
    [dates, discounts1] = bootstrap(datesSet, shiftedRates);
    discounts1_reset = interpolation(discounts1, dates, resetDates_all);
    discounts1_reset(1) = 1;

    fwd_discounts = discounts1_reset(2:end)./discounts1_reset(1:end-1);
    deltas = yearfrac(resetDates_all(1:end-1),resetDates_all(2:end),ACT_360);

    euribor1 = 1./deltas.*(1./fwd_discounts-1);
    euribor1 = euribor1(2:end);
   
    sigma_new = calibrateLMM(euribor1, mkt_strikes, resetDates_all, mkt_vols, dates, discounts1);
    sigma_new =  sigma_new(1:40,:);
    euribor1 = euribor1(1:40);

    upfront1 = computeUpfront(euribor1, mkt_strikes, sigma_new, spol, resetDates, dates, discounts1, "Black");

    upfrontDiff(i)=upfront1-upfront_black;
end
toc;

fprintf("Delta-Sens 2y: %d\n", upfrontDiff(13));
fprintf("Delta-Sens 6y: %d\n", upfrontDiff(17));
fprintf("Delta-Sens 10y: %d\n", upfrontDiff(end));

%% Point 1.4 - Total Vega

% We shift upwards the volatility surface by 1bp and compute the difference
% in terms of upfront. This is the contract sensitivity to the Vega.

bp=1e-4;
mkt_vols_shocked = mkt_vols + bp;

discountsReset_all = interpolation(discounts, dates, resetDates_all);
discountsReset_all(1) = 1;
fwd_discounts_all = discountsReset_all(2:end)./discountsReset_all(1:end-1);
deltas_all = yearfrac(resetDates_all(1:end-1),resetDates_all(2:end),ACT_360);

euribor_all= 1./deltas_all.*(1./fwd_discounts_all-1);
euribor_all = euribor_all(2:end);
sigma_shocked = calibrateLMM(euribor_all, mkt_strikes, resetDates_all, mkt_vols_shocked, dates, discounts);

upfront_shocked = computeUpfront(euribor, mkt_strikes,sigma_shocked(1:40,:),spol,resetDates,dates,discounts,"Black");
upfrontDiffSigma = upfront_shocked - upfront_black;

fprintf("Vega (Single contract): %d\n", upfrontDiffSigma);

%% Point 1.5 - Course Grained Delta Hedging (2y, 6y, 10y)

ratesSet_buck2y = bucket_2y(ratesSet);
[~, discounts2y] = bootstrap(datesSet, ratesSet_buck2y);
sigma_new = calibrateLMM(euribor_all, mkt_strikes, resetDates_all, mkt_vols, dates, discounts2y);
upfront2y= computeUpfront(euribor,mkt_strikes,sigma_new(1:40,:),spol,resetDates,dates,discounts2y,"Black");
delta_2y = upfront2y - upfront_black;
ratesSet_buck6y = bucket_6y(ratesSet);
[~, discounts6y] = bootstrap(datesSet, ratesSet_buck6y);
sigma_new = calibrateLMM(euribor_all, mkt_strikes, resetDates_all, mkt_vols, dates, discounts6y);
upfront6y= computeUpfront(euribor,mkt_strikes,sigma_new(1:40,:),spol,resetDates,dates,discounts6y,"Black");
delta_6y = upfront6y - upfront_black;
ratesSet_buck10y = bucket_10y(ratesSet);
[~, discounts10y] = bootstrap(datesSet, ratesSet_buck10y);
sigma_new = calibrateLMM(euribor_all, mkt_strikes, resetDates_all, mkt_vols, dates, discounts10y);
upfront10y= computeUpfront(euribor,mkt_strikes,sigma_new(1:40,:),spol,resetDates,dates,discounts10y,"Black");
delta_10y = upfront10y - upfront_black;

swapRate2y= computeSwapRate(dates, discounts, resetDates(1:4:9));
swapRate6y= computeSwapRate(dates, discounts, resetDates(1:4:25));
swapRate10y= computeSwapRate(dates, discounts, resetDates(1:4:end));

ratesSetShifted = ratesShiftAll(ratesSet);
[~, discountsShifted] = bootstrap(datesSet, ratesSetShifted(1:4:end));
sensSwap2y = NPV_irs(dates, discountsShifted, resetDates(1:4:9) ,swapRate2y);
sensSwap6y = NPV_irs(dates, discountsShifted, resetDates(1:4:25) ,swapRate6y);
sensSwap10y = NPV_irs(dates, discountsShifted, resetDates(1:4:end),swapRate10y);

notional = 50e6;

[x2,x6,x10] = hedgeCertificate(sensSwap2y, sensSwap6y, sensSwap10y, delta_2y, delta_6y, delta_10y, notional);

fprintf("--------Delta HEDGING--------\n");
fprintf("IRS 2y: %d\n", x2);
fprintf("IRS 6y: %d\n", x6);
fprintf("IRS 10y: %d\n", x10);

%% Point 1.6 - Course Grained Vega Hedging (6y, 10y)

mkt_vols_6y = vegaBucket_6y(mkt_vols);
mkt_vols_10y = vegaBucket_10y(mkt_vols);

sigmaShocked6y = calibrateLMM(euribor_all, mkt_strikes, resetDates_all, mkt_vols_6y, dates, discounts);
sigmaShocked10y = calibrateLMM(euribor_all, mkt_strikes, resetDates_all, mkt_vols_10y, dates, discounts);

vega_6y = computeUpfront(euribor, mkt_strikes, sigmaShocked6y(1:40,:), spol, resetDates, dates, discounts,"Black") - upfront_black;
vega_10y = computeUpfront(euribor, mkt_strikes, sigmaShocked10y(1:40,:), spol, resetDates, dates, discounts,"Black") - upfront_black;

capVega_6y = vegaCap(euribor, mkt_strikes, resetDates(1:25), sigma, dates, discounts, swapRate6y);
capVega_10y = vegaCap(euribor, mkt_strikes, resetDates, sigma, dates, discounts, swapRate10y);

[v6,v10] = hedgeVega(vega_6y, vega_10y, capVega_6y, capVega_10y, notional);

fprintf("--------Vega HEDGING--------\n");
fprintf("CAP 6y: %d\n", v6);
fprintf("CAP 10y: %d\n", v10);

%% Point 2.1 - Calibrating the Bond Market Model (BMM)

tic;
sigma_BMM = calibrateBMM(euribor, mkt_strikes, resetDates, mkt_vols, dates, discounts);
toc;

[T,K] = meshgrid(resetDates(2:end-1), mkt_strikes);

surf(T,K,sigma_BMM')
title("Bond Market Model - Volatility Surface");

%% Point 2.2 - Pricing the option

lambda=0.1;
B_matrix = B_MC(dates, discounts, resetDates(1:17), sigma_BMM, lambda, euribor, mkt_strikes);
[n,m]=size(B_matrix);
L_matrix=B_matrix;
ACT_365=3;
yf=yearfrac(resetDates(1:16), resetDates(2:17), ACT_365);

for i=1:m
    L_matrix(:,i) = 1/yf(i) * (1./B_matrix(:,i) - 1);
end

payoff=[];

for i=1:length(resetDates(3:17))
    payoff = [payoff, yf(i+1) * max(0, L_matrix(:,i+1)-L_matrix(:,i)-0.0005)];
end

%% Point 2.2 - Pricing the option

discounts_inter = interpolation(discounts,dates,resetDates(3:17));
priceBMM=zeros(1e7,1);
for i=1:1e7
    priceBMM(i)=sum(discounts_inter.*payoff(i,:));
end
stDev=std(priceBMM);
confLevel=0.99;
priceBMM =mean(priceBMM);
CI =[priceBMM - norminv(1-confLevel/2)*stDev, priceBMM + norminv(1-confLevel/2)*stDev];

fprintf("Price of the option: %d\n", priceBMM);
fprintf("Confidence Interval: [%0.4f,%0.4f]\n", CI(1), CI(2));