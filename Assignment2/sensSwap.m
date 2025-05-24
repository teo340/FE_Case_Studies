function [DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts,discounts_DV01)
% Sensitivity of the Interest Rates Swap.
% This function computes the DV01, BPV and DV01_z of the swap. 
% INPUT
% setDate: settlement date of the swap, format datetime
% fixedLegPaymentDates: payment dates of the fixed leg, format datetime
% fixedRate: swap rate
% dates: dates of the curve. We ended up not using this parameter
% discounts: discount factor curve
% discounts_shifted: discount factor curve shifted by 1bp
% OUTPUT
% DV01: Dollar Value of 1 basis point
% BPV: Basis Point Value
% DV01_z: Dollar Value of 1 basis point, zero rates

% We set t0 as the settlement date
t0 = setDate;
EU_30_360 = 6; % daycount convention

% Forward deltas - year fraction between payment dates
yf = yearfrac(t0, fixedLegPaymentDates, EU_30_360);
deltas = yearfrac(t0, fixedLegPaymentDates(1), EU_30_360);
deltas = [deltas; yf(2:end)-yf(1:end-1)];

% As already mentioned in the bootstrap function, the discount factor
% at 1 year is not included in the curve. However, it is needed right now.
% We chose to hard code it here, to not modify the function signature provided
% by the assignment. To see how this number was computed, check out
% the bootstrapSwap function.
discounts = [0.968301310232409; discounts(13:18)]; % discounts for original curve
discounts_DV01 = [0.968203894341989; discounts_DV01(13:18)]; % discounts for shifted curve

bp = 0.0001; % 1 basis point

% BPV
BPV = sum(deltas.*discounts)*bp;

% DV01
NPV_0 = 1-discounts(end) - fixedRate*sum(deltas.*discounts); % pratically zero - in accordance to the theory
NPV_shifted = 1-discounts_DV01(end) - fixedRate*sum(deltas.*discounts_DV01);
DV01 = abs(NPV_shifted-NPV_0);

% DV01_z 
ACT_365 = 3; % daycount convention
zeroRates_0 = zeroRates(fixedLegPaymentDates, discounts);
zeroRates_shifted = zeroRates_0+bp;
discounts_shifted = exp(-zeroRates_shifted.*yearfrac(t0, fixedLegPaymentDates, ACT_365));

NPV_z = 1 - discounts_shifted(end) - fixedRate * sum(deltas.*discounts_shifted);
DV01_z = abs(NPV_z-NPV_0);

end     % function sensSwap