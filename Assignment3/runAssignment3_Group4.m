% runAssignment3
% group 4, AY2024-2025
% Based on the bootstrap computed in the previous assignment, we
% compute some quantities like the asset swap rate or the cds spreads 
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

%% Point 1
% We consider an asset swap with the following characteristics:
% - Coupon: 4.8%
% - start date: 02-Feb-2022
% - end date: 02-Feb-2028
% - today: 02-Feb-2023
% - clean price: 101
% We compute the asset swap rate

coupon = 0.048;      % Coupon of the corresponding bond
% We generate a set of dates for the fixed payment dates of the bond
startYear = 2022; 
endYear = 2028;
today = datetime('02-Feb-2023');
cleanPrice = 101;    % Clean price of the bond, no accrual

% Computation of the asset swap rate with the function 'asset_swap'
as_rate = asset_swap(discounts, dates, coupon, cleanPrice, startYear, endYear, today);
fprintf('Asset Swap Rate: %d\n', as_rate);

%% Point 2
% Define parameters
bp = 1e-4;
spreadsCDS = [40, 44, 47, 49, 51, 52]*bp;
datesCDS = datetime(2024:2029,2,2);
flag = 1;
recovery = 0.4;

% check business days
datesCDS = check_busday(datesCDS);

[datesCDS, survival, intesities] = bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, 1, recovery);
[~, survival_acc, intesities_acc] = bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, 2, recovery);
[~, survival_jt, intesities_jt] = bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, 3, recovery);

% Plot of the survival probabilities 
figure;
plot(datesCDS, survival(2:end));
hold on; grid on;
plot(datesCDS, survival_acc(2:end),'--');
plot(datesCDS, survival_jt,'-*');
legend('Survival','with accrual','JT');

hold off
figure
hold on; grid on
% Plot of the intensities piecewise constant 
x_int = [datesCDS(1:end-1); datesCDS(2:end)];
y_int = [intesities(1:end-1)'; intesities(1:end-1)'];
x_int = [x_int, [datesCDS(end); datesCDS(end) + (datesCDS(end) - datesCDS(end-1))]];  
y_int = [y_int, [intesities(end); intesities(end)]];
plot(x_int, y_int, 'b', 'LineWidth', 2);

title('Intensities');
legend('Intensities (neglecting the accrual)');


%% Point 3
% Define parameters
bp = 1e-4;
lambdas = [5,9]*bp;
theta = 4; % years
M = 1e5;
% Compute a sample of default times with the function 'credit_simulation'
tau_samples = credit_simulation(lambdas, theta, M);
% Compute an estimation of the parameters and the Confidence Interval
[lambdas_est, CI] = estimate_params(tau_samples, theta);

% Print the results
fprintf('lambda_1: %d\n', lambdas_est(1));
fprintf('lambda_2: %d\n', lambdas_est(2));
disp('-------------')

%disp(table(beta, CI_lower, CI_upper, 'VariableNames', {'Beta', 'LowerCI', 'UpperCI'}));
disp(CI)


%% Point 4
% Define parameters
faceValue = 1e9;
correlation = 0.4;
recovery = 0.2;
p = 0.05;
% we had to run the program with a limited number of bond holders, 
% because the exact solution presented numerical problems for larger I. 
I = 100; 
Kd = 0.05; Ku = 0.09;
N_mezz=faceValue*(Ku-Kd);

% Pricing of the MBS
price = mbs_pricing(dates, discounts, I, p, correlation, Ku, Kd, recovery, 1);
price_2 = mbs_pricing(dates, discounts, I, p, correlation, Ku, Kd, recovery, 2);
price_3 = mbs_pricing(dates, discounts, I, p, correlation, Ku, Kd, recovery, 3);

% Print the results
fprintf('HP: %d (%.2f%%) \n', price*N_mezz, price*100);
fprintf('KL: %d (%.2f%%)\n', price_2*N_mezz, price*100);
fprintf('LHP: %d (%.2f%%)\n', price_3*N_mezz, price*100);

%% This code snippet may take around 10 seconds to run!

tic
N=10;
I = logspace(1,5,N);
prices = zeros(2,N);
prices_HP = zeros(1,4);
for i=1:N
    if i <= 4
        prices_HP(i) = mbs_pricing(dates, discounts, floor(I(i)), p, correlation, Ku, Kd, recovery, 2);
    end

    prices(1,i) = mbs_pricing(dates, discounts, floor(I(i)), p, correlation, Ku, Kd, recovery, 2);
    prices(2,i) = mbs_pricing(dates, discounts, i, p, correlation, Ku, Kd, recovery, 3);
end
toc

figure;
semilogx(floor(I(1:4)), prices_HP*N_mezz, 'b-o', 'LineWidth', 4);
hold on; grid on;
semilogx(floor(I), prices(1,:)*N_mezz, 'r-s', 'LineWidth', 2);
semilogx(floor(I), prices(2,:)*N_mezz, 'g-d', 'LineWidth', 2);
legend('HP', 'KL', 'LHP');
hold off

%% Point b - HP and LHP

Kd = 0.0; Ku = 0.05;
I = logspace(1,5,10);
N_equity = (Ku-Kd)*faceValue; 
prices = zeros(2,N);
prices_HP = ones(1,4);
for i=1:N
    if i <= 4
        [prices_HP(i),~] = mbs_pricing(dates, discounts, floor(I(i)), p, correlation, Ku, Kd, recovery, 1);
    end
    [prices(2,i),~] = mbs_pricing(dates, discounts, I(i), p, correlation, Ku, Kd, recovery, 3);
end

figure;
semilogx(floor(I(1:4)), prices_HP*N_equity, 'b-o', 'LineWidth', 2);
hold on; grid on;
semilogx(floor(I), prices(2,:)*N_equity, 'r-s', 'LineWidth', 2);
legend('HP','LHP');

%% Point B - KL

Kd_ME = 0.0; Ku_ME = 0.09; % mezzanine + equity
[price_ME,~] = mbs_pricing(dates, discounts, 100, p, correlation, Ku_ME, Kd_ME, recovery, 2);
notional_ME = (Ku_ME-Kd_ME)*1e9;
price_ME = notional_ME*price_ME;

Kd_M = 0.05; Ku_M = 0.09;
[price_M, ~] = mbs_pricing(dates, discounts, 100, p, correlation, Ku_M, Kd_M, recovery, 2);
notional_M = (Ku_M-Kd_M)*1e9;
price_M = price_M*notional_M;

Ku_E = 0.05; Kd_E = 0;
notional_E = (Ku_E-Kd_E)*1e9;

price_E = price_ME - price_M;
fprintf('Price of the Equity Tranch: %d\n', price_E);

