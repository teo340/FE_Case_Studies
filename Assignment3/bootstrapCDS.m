function [datesCDS, survProbs, intensities] = bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)
% This function computes the survival probabilities and intensities implied by the CDS spreads.
% We assume the intensities are piece-wise constant between the CDS dates.
% INPUTS:
% datesDF: dates of the discount factors, bootstrapped
% discounts: discount factors, bootstrapped
% datesCDS: dates of the CDS contracts, format datetime
% spreadsCDS: CDS spreads, in basis points
% flag: 1 -> no accrual, 2 -> with accrual, 3 -> Jarrow-Turnbull approximation
% recovery: recovery rate
% OUTPUTS:
% datesCDS: dates of the CDS contracts, format datetime
% survProbs: survival probabilities
% intensities: intensities

% We define the convention we will use in the code (30/360 European convention)
EU_30_360 = 6;

% Number of CDS maturities + 1 for today    
N = length(datesCDS)+1;
% Initialize survival probabilities (P) with ones
P = ones(N,1);

% Define the reference date (today)
today = datetime('02-Feb-2023');
% Interpolate discount factors at CDS dates
discounts = interpolation(discounts, datesDF, datesCDS);
% We compute time fractions between key dates
deltas = [
    yearfrac(today, datesCDS(1), EU_30_360);
    yearfrac(datesCDS(1:end-1), datesCDS(2:end), EU_30_360)';
];

% Case 1: No accrual
if flag == 1
    for i=2:N
     % We compute the first sum term in the survival probability formula
        sum1 = (1-recovery)*sum(discounts(1:i-1).*P(1:i-1)');

     % We compute the second sum term (only if i > 2)
        if i>2
            sum2 = sum(discounts(1:i-2).*P(2:i-1)'.*(spreadsCDS(i-1)*deltas(1:i-2)'+(1-recovery)));
        else
            sum2 = 0;
        end

        % We evaluate survival probability P(i) according to the theory
        num = sum1-sum2;
        den = discounts(i-1)*(spreadsCDS(i-1)*deltas(i-1)+(1-recovery));
    
        P(i) = num/den;
    end

    % We store survival probabilities
    survProbs = P;
    % Compute intensities using actual/365 day count convention
    ACT_365 = 3;
    deltas = [
        yearfrac(today, datesCDS(1), ACT_365);
        yearfrac(datesCDS(1:end-1), datesCDS(2:end), ACT_365)';
    ];
    intensities = -1./deltas.*log(P(2:end)./P(1:end-1));

% Case 2: With accrual
elseif flag==2

    for i=2:N
        % We compute the first sum term with accrual adjustment
        sum1 = sum(discounts(1:i-1).*P(1:i-1)'.*((1-recovery)-0.5*spreadsCDS(i-1)*deltas(1:i-1)'));
    
        if i>2
            % We calculate the second sum term (only if i > 2 as before)
            sum2 = sum(discounts(1:i-2).*P(2:i-1)'.*(0.5*spreadsCDS(i-1)*deltas(1:i-2)'+(1-recovery)));
        else
            sum2 = 0;
        end

        % Compute survival probability P(i)
        num = sum1-sum2;
        den = discounts(i-1)*(0.5*spreadsCDS(i-1)*deltas(i-1)+(1-recovery));
        P(i) = num/den;
    end

    % Store survival probabilities
    survProbs = P;
    % We compute intensities using actual/365 day count convention
    ACT_365 = 3;
    deltas = [
        yearfrac(today, datesCDS(1), ACT_365);
        yearfrac(datesCDS(1:end-1), datesCDS(2:end), ACT_365)';
    ];
    intensities = -1./deltas.*log(P(2:end)./P(1:end-1));

% Case 3: Jarrow-Turnbull approximation
elseif flag==3
    % Jarrow-Turnbull approximation
    ACT_365 = 3;
    
    % We compute year fractions from today to CDS dates
    yf = yearfrac(today, datesCDS, ACT_365);
    
    % We evaluate intensities directly using spreads and recovery rate
    intensities = spreadsCDS./(1-recovery);
   
    % We calculate survival probabilities according to the theory
    survProbs = exp(-intensities.*yf);
end

end
