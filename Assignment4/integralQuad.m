function I = integralQuad(f,logMoneyness,N)
% This function numerically integrates the function f,
% over the space range -N to N, using the quadrature method.
% INPUTS:
% f: function to be integrand (without exp(-i*u*z)
% logMoneyness: vector of logMoneyness values
% N: integration extrema
% OUTPUT:
% I: vector of integral values at the logMoneyness points

I = zeros(length(logMoneyness),1);

for j=1:length(logMoneyness)
    integrand = @(u) exp(-1i*u.*logMoneyness(j)).*f(u);
    I(j) = quadgk(integrand,-N,N);
end

end