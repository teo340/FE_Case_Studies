function [price] = callPricing(logMoneyness, sigma, k, eta, alpha, F0, B, t, flag)
% This function prices a call option for every logMoneyness,
% using a wide variety of methods: FFT, Quadrature or Monte Carlo. The
% returns of the forward are assumed to follow a Normal Mean Variance Mixture model.
% INPUTS:
% logMoneyness: vector of logMoneyness values
% sigma: level of volatility
% k: volatility of volatility
% eta: skewness
% alpha: NMVM parameter
% F0: forward price
% B: discount factor
% t: maturity date
% flag: method to use for pricing ('FFT', 'Quad', 'MC')
% OUTPUTS:
% price: vector of call option prices for each logMoneyness value

% Depending on alpha, the laplace transform of the characteristic function
% changes. Here we define two cases: alpha = 0 and alpha > 0.
if alpha == 0
    % G is a Gamma distribution with parameters a and b
    laplaceExponent = @(u) -t/k*log(1+k*u*sigma^2);
else
    % G is an Invese Gaussian
    laplaceExponent = @(u) t/k*(1-alpha)/alpha*(1-(1+(u*k*sigma^2)./(1-alpha)).^alpha);
end
   
% Characteristic function of the forward price process
mu = -laplaceExponent(eta);
phi = @(u) exp(1i*u*mu).*exp(laplaceExponent((u.^2+1i*(2*eta+1)*u)/2));

% FFT and Quad parameters
M = 15;
dz = 0.001;

if flag == "Quad"
    f = @(u) phi(-u-1i/2)./(u.^2+1/4);
    I = integralQuad(f, logMoneyness, 10000);
    I = real(I);
    price = 1-exp(-logMoneyness'/2)/(2*pi).*I;
    price = F0*B*price;
elseif flag == "FFT"
    I = integralFFT(phi, M, dz, logMoneyness);
    price = 1-exp(-logMoneyness/2).*I;
    price = F0*B*price;
elseif flag == "MC"
    [price, IC] = integralMC(sigma, k, eta, t, logMoneyness, F0, B);
   
    mask = (logMoneyness == 0);
    fprintf("Confidence Interval ATM: [%.2f, %.2f]\n\n", IC(1, mask), IC(2, mask));

    figure;
    plot(logMoneyness, IC(1,:),'red');
    hold on; grid on;
    plot(logMoneyness, IC(2,:),'red');
    plot(logMoneyness, price, 'o');
    hold off

end

end