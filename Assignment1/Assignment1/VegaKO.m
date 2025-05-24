function vega = VegaKO(F0,K,KO,B,T,sigma,N,flagNum)
% This function computes the Vega of a Knock-Out option,
% using either the closed solution, the CRR approach or
% the Monte Carlo Method. 
% INPUTS
% F0: Initial price of the underlying asset
% K: Strike price
% KO: Knock-Out barrier
% B: Discount factor
% T: Maturity
% sigma: Volatility
% N: Number of steps for CRR, number of simulations for MC
% flagNum: 1 for CRR, 2 for MC, 3 for closed form
% OUTPUT
% vega: Vega of the Knock-Out option

switch flagNum
    case 1
        % CRR
        vega = VegaKOCRR(F0,K,KO,B,T,sigma,N);
    case 2
        % MC
        vega = VegaKOMC(F0,K,KO,B,T,sigma,N);
    case 3
        % We implement the Vega in closed form. Because of linearity,
        % the Vega of the Knock-Out option is the sum of the vega of the
        % first call, minus the vega of the second, minus the vega of the
        % digital option. The latter vega has been calculated by hand.

        d1_CallK = log(F0/K)/(sigma*sqrt(T))+0.5*sigma*sqrt(T);
        d1_CallKO = log(F0/KO)/(sigma*sqrt(T))+0.5*sigma*sqrt(T);

        d2_DigKO = log(F0/KO)/(sigma*sqrt(T))-0.5*sigma*sqrt(T);
        d1_DigKO = log(F0/KO)/(sigma*sqrt(T))+0.5*sigma*sqrt(T);

        % Call K
        vega1 = F0*B*sqrt(T)*exp(-d1_CallK^2/2)/sqrt(2*pi); % vega of the first call
        vega2 = F0*B*sqrt(T)*exp(-d1_CallKO^2/2)/sqrt(2*pi); % vega of the seconda call
        vega3 = -B/sigma*d1_DigKO*normpdf(d2_DigKO); % vega of the digital option

        dSigma = 0.01; % 1% of basis point
        vega = (vega1-vega2-(KO-K)*vega3)*dSigma;
    otherwise
        error('Invalid flagNum. Choose 1 for CRR, 2 for MC, or 3 for closed form.');
end

end     % function VegaKO