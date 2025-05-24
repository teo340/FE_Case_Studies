function b = bin_LHP(I,m)
% This function computes the binomial coefficient
% using the Stirling-De Moivre approximation.
% For large numbers, it doesn't overflow.
% INPUTS:
% I: number of trials
% m: number of successes
% OUTPUTS:
% b: binomial coefficient

    if m == I
        b = 1;
        return;
    end
% 
%     z = m/I;
%     b = 1./sqrt(I).*1./sqrt(2*pi*z.*(1-z)).*z.^(-z.*I).*(1-z).^(-I.*(1-z));
    z = m./I;
    C = sqrt(I./(2*pi*(1-z).*z));
    b = 1./I.*C.*exp(-I.*(z.*log(z)+(1-z).*log(1-z)));
end     % function bin_LHP