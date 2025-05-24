function [price,expected_loss] = mbs_pricing(dates, discounts, I, p, correlation, Ku, Kd, recovery, flag)
% The function computes the price and expected loss of a mortgage-backed security (MBS)
% using different pricing methodologies based on the specified flag for a mezzanine tranche.
%
% Inputs:
%   dates: vector of dates corresponding to discount factors.
%   discounts: vector of discount factors.
%   I: number of obligors.
%   p: probability of default.
%   correlation: correlation factor for default dependency.
%   Ku: upper detachment point in percentage.
%   Kd: lower detachment point in percentage.
%   recovery: recovery rate.
%   flag: method selection: 1 (HP), 2 (KL), 3 (LHP).
%
% Outputs:
%   price: computed MBS price.
%   expected_loss: expected loss computed based on the selected method.

% We define the needed parameters for our model
d = Kd/(1-recovery);
u = Ku/(1-recovery);

% We define the loss function according to the theory
loss = @(z) min(max(z-d,0),u-d)/(u-d); 
K = norminv(p);

% We define the obligor default probability for a given y
p_y = @(y) normcdf((K-sqrt(correlation)*y)./sqrt(1-correlation));

% Initialize expected loss
expected_loss = 0;

% First case: Homogeneous portfolio
if flag==1
    for m=0:I
        % We compute the expected loss with a loop circle 
        expected_loss = expected_loss + loss(m/I)*prob_default(I, m, p, correlation, 1);
    end
% Second case: Kullback-Leibler
elseif flag==2
    % We compute the needed quantities, included the Kullback-Leibler entropy (function of two variables: z, y)
    K = norminv(p);
    p_y = @(y) normcdf((K-sqrt(correlation)*y)./sqrt(1-correlation));
    K_entropy = @(z,y) z.*log(z./p_y(y))+(1-z).*log((1-z)./(1-p_y(y)));

    % Normalization function of the dentity of z given y derived from Stirling formula    
    C1 = @(z) sqrt(I./(2*pi*z.*(1-z)));

    % We define the integrand function
    integrand = @(z,y) C1(z).*exp(-I*K_entropy(z,y));
    % The matlab command 'arrayfun' is necessary in order to make the code work since we are dealing with vectors as input, 
    % in this way the handle function D is capable of accepting vectors
    D = @(y) arrayfun(@(yy) integral(@(z) integrand(z, yy), 0, 1), y);
    C = @(z,y) C1(z)./D(y);

    % We compute the integrals for expected loss
    integrand1 = @(z,y) loss(z).*C(z,y).*exp(-I.*K_entropy(z,y));
    quad_integral = @(y) arrayfun(@(yy) quadgk(@(z) integrand1(z, yy), 0, 1), y);
    integrand2 = @(y) normpdf(y).*quad_integral(y);
    expected_loss = integral(integrand2,-6,6);

% Third case: Large Homogeneous Portfolio
elseif flag==3

    K = norminv(p);
    p_y = @(y) normcdf((K-sqrt(correlation)*y)./sqrt(1-correlation));
    
    % We define the inverse of the function p
    inv_p = @(z) (K - sqrt(1-correlation).*norminv(z))./sqrt(correlation);
    
    % We compute the first derivative of the previous function
    der_inv_p = @(z) 1./normpdf(norminv(z)).*(sqrt((1-correlation)/correlation));
    
    % We define the density function for the fraction of defaulted obligors z=m/I 
    f_LHP= @(z) normpdf(inv_p(z)).*der_inv_p(z);
    
    % Compute expected loss using integral
    integrand=@(z) loss(z).*f_LHP(z);
    expected_loss=integral(integrand,0,1);

end

% We define the maturity
maturity = datetime('02-Feb-2026');
% We compute the required discount factor to evaluate the price 
B = interpolation(discounts, dates, maturity);

% We Compute the final price
price = B*(1-expected_loss);


end
