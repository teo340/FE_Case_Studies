% Assignment_1
%  Group 4, AA2024-2025
%
clear;
close;
clc;

%% Pricing parameters
S0=1;
K=1.05;
r=0.025;
TTM=1/3; 
sigma=0.21;
flag=1; % flag:  1 call, -1 put

d=0.02;
NumberOfContract = 1e6;

%% Quantity of interest
B=exp(-r*TTM); % Discount

%% Pricing 
F0=S0*exp(-d*TTM)/B;     % Forward in G&C Model

rng("default"); % set random seed for MC reproducibility

% We price the option according to three different methods:
% 1 ClosedFormula, 2 CRR, 3 Monte Carlo
prices = zeros(3,1);
M=100; % M = simulations for MC, steps for CRR;
for pricingMode=1:3
    prices(pricingMode) = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag);
    fprintf('Option price: %0.8f\n', prices(pricingMode));
    fprintf('Option price: %0.8f\n', 1e6*prices(pricingMode));
end
fprintf('-------------\n');

%% Errors Rescaling 
% plot Errors for CRR varing number of steps
% Note: both functions plot also the Errors of interest as side-effect 
[nCRR,errCRR]=PlotErrorCRR(F0,K,B,TTM,sigma);

% plot Errors for MC varing number of simulations N 
[nMC,stdEstim]=PlotErrorMC(F0,K,B,TTM,sigma); 

spread = 1e-4;
for i=1:length(errCRR)
    if errCRR(i) <= spread
        fprintf('CRR converged in %d steps\n',nCRR(i));
        break;
    end
end
% with this seed, the variance of the first 3 MC simulations
% is very small, while the error in the price is still high.
% We decided to not consider the first 3.
for i=3:length(stdEstim)
    if stdEstim(i) <= spread
        fprintf('MC converged in %d simulations\n',nMC(i));
        break;
    end
end

M_optimal = findM(prices(1), F0, K, B, TTM, sigma, 1, 1e-4);
fprintf('Optimal numnber iterations: %d\n', M_optimal);


%% KO Option - Barrier Option
% replicated by LongCall(K) Short Call(KO) Short Digital(KO-K)
% We price a KnockOut option with a KO = 1.4
KO = 1.4;
M = 100;
KOprice_closed = EuropeanOptionKO(F0,K,KO,B,TTM,sigma);
KOprice_mc = EuropeanOptionKOMC(F0,K,KO,B,TTM,sigma,1000);
KOprice_crr = EuropeanOptionKOCRR(F0,K,KO,B,TTM,sigma,M);

fprintf('KO price - Closed solution: %0.6f\n', KOprice_closed);
fprintf('KO price - Monte Carlo:     %0.6f\n', KOprice_mc);
fprintf('KO price - CRR:             %0.6f\n', KOprice_crr);
fprintf('-------------\n');

%% KO Option Vega
% We plot the option Vega in three different ways,
% according to the closed form solution, the CRR approximation
% and the MC simulation. Derivatives have been approximated
% using a centered finite difference scheme.
% The function may take some time to run, for large
% numbers of iterations. Here we use M = 5000

price_space = 0.65:0.01:1.45;
vega_closed = zeros(1,length(price_space));
vega_crr = zeros(1,length(price_space));
vega_mc = zeros(1,length(price_space));
for i=1:length(price_space)
    F_i = price_space(i)*B;
    vega_crr(i) = VegaKO(F_i,K,KO,B,TTM,sigma,5000,1);
    vega_mc(i) = VegaKO(F_i,K,KO,B,TTM,sigma,5000,2);
    vega_closed(i) = VegaKO(F_i,K,KO,B,TTM,sigma,5000,3);
end

figure
subplot(1,2,1);
plot(price_space, vega_closed);
hold on; grid on;
plot(price_space, vega_crr);
%ylim([-5e-3,3e-3]);
legend('Closed','CRR');
title('Vega: Closed vs CRR')

subplot(1,2,2);
plot(price_space, vega_closed);
grid on; hold on;
plot(price_space, vega_mc);
%ylim([-5e-3,5e-3]);
legend('Closed','MC');
title('Vega: Closed vs MC')

%% Antithetic Variables
% The use of antithetic variables can significantly reduce the number
% of MonteCarlo iterations need to converge.
% We compute the standard deviation of prices for both: the standard
% MonteCarlo and the Varianced Reduced(_vr) version.
m = 1:20;
M = 2.^m;
stdEstim = zeros(length(M),1);
stdEstim_vr = zeros(length(M),1);

% we generate M(i) samples of a standard gaussian random variable, 
% to then compute the discounted payoff.
for i=1:length(M)
    prices_1 = B*max(F0*exp(-sigma^2/2*TTM+sigma*sqrt(TTM)*rand(M(i),1))-K,0);
    prices_2 = B*max(F0*exp(-sigma^2/2*TTM-sigma*sqrt(TTM)*rand(M(i),1))-K,0);
    prices_vr = (prices_1+prices_2)/2; % antithetic vars technique

    prices = B*max(F0*exp(-sigma^2/2*TTM+sigma*sqrt(TTM)*rand(M(i),1))-K,0);

    stdEstim(i) = 1/sqrt(M(i))*std(prices);
    stdEstim_vr(i) = 1/sqrt(M(i))*std(prices_vr);
end

figure;
loglog(M,stdEstim,'LineWidth',1.5);
hold on; grid on;
loglog(M,stdEstim_vr,'LineWidth',1.5);
loglog(M,1./sqrt(M),'LineWidth',1.5);
loglog(M,spread*ones(length(M),1),'LineWidth',1.5,'LineStyle','--');
title('Error MC');
legend('MC','MC-AV','1/sqrt(M)','spread');
hold off;
% we see from the plot that the spread is reached
% in almost half order of magnitude. 

%% Bermudan Option: Pricing
% We price a Bermudan option, with exercise possibilities
% at the end of each month. We do so using a Tree based approach.
% We set the number of time steps N. For now, we quote with a dividend
% rate d = 0.02

M = 100;
ber_price = BermudanOptionCRR(F0, K, B, TTM, sigma, d, M);
fprintf('Bermudan price: %0.8f\n', ber_price);

% The Bermudan and the European are very close, because the dividend yield
% is still small. As we'll see shortly, the price difference will widen
% as we increase the dividend.

%% Bermudan Option: Dividend Yield
% We study how the price how the Bermudan deviates from the price
% of a Vanilla european. We consider dividends in the range (0,5%).
d_step = 100;
d_space = linspace(0,0.05,d_step);

% We compute the price of the European and of the Bermudan for each value
% of the dividend and plot the result.
eu_price = zeros(1,d_step);
ber_price = zeros(1,d_step);
for i=1:d_step
    F_i =S0*exp(-d_space(i)*TTM)/B;     % Forward in G&C Model
    eu_price(i) = EuropeanOptionCRR(F_i,K,B,TTM,sigma,M,1);
    ber_price(i) = BermudanOptionCRR(F_i,K,B,TTM,sigma,d_space(i),M);
end

figure;
plot(d_space, eu_price,'LineWidth',1.2,'LineStyle',':');
hold on; grid on;
plot(d_space, ber_price,'LineWidth',1.5,'LineStyle','--');
legend('European','Bermudan');
xlabel('Dividend Yield');
ylim([0.025 0.03]);
title('Dividend vs Bermudan');
