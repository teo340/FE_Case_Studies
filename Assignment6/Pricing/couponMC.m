function [coupon,Ft] = couponMC(F0, strike, sigma, k, eta, today, T, B)
% The function simulates the underlying trajectory, according
% to a Variance-Gamma process and computes the coupon of the
% certificate.
% INPUTS:
% F0: initial forward price
% strike: strike of the certificate
% sigma: volatility of the forward price process
% k: parameter of the Variance-Gamma process
% eta: parameter of the Variance-Gamma process
% today: today's date
% T: maturity date of the certificate
% B: vector of discount factors for the payment dates

ACT_365 = 3; 
N_sim = 1e5;
rng(42);
resetDates = businessdayoffset(today:calyears(1):T);
deltaT = yearfrac(resetDates(1:end-1), resetDates(2:end), ACT_365);

Fprev = F0;
Ft = zeros(N_sim, 3);
Ft(:,1) = F0*ones(N_sim,1);
for i=1:length(resetDates)-1
    g = randn(N_sim,1);
    a = deltaT(i)/k; b = 1/a;
    G = gamrnd(a,b,[N_sim,1]);

    % Characteristic function of the forward price process
    laplaceExponent = @(u) -deltaT(i)/k*log(1+k*u*sigma^2);
    mu = -laplaceExponent(eta);
    f = mu - (0.5+eta)*sigma^2*G*deltaT(i)+sigma*sqrt(deltaT(i))*sqrt(G).*g;
    
    Ft(:,i+1) = Fprev.*exp(f);
    
    Fprev = Ft(:,i+1);
end

EU_30_360 = 6;
delta_1y = yearfrac(today, resetDates(1), EU_30_360);
delta_2y = yearfrac(resetDates(1), resetDates(2), EU_30_360);
delta_3y = yearfrac(resetDates(2), resetDates(3), EU_30_360);

% The coupon is not fixed, but depends on the trajectory of the underlying.
% The first coupon is always paid, while the second coupon is paid
% only if the early redemption clause is not triggered. Same goes for the third coupon.

coupon = B(1)*delta_1y*0.06*(Ft(:,2) < strike) + ...
         B(2)*delta_2y*0.06*(Ft(:,2) >= strike & Ft(:,3) < strike) + ...
         B(3)*delta_3y*0.02*(Ft(:,2) >= strike & Ft(:,3) >= strike);

end