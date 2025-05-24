function as_rate = asset_swap(discounts, dates, coupon, cleanPrice, startYear, endYear, today)
% This function computes the asset swap rate
% INPUTS:
% discounts: discounts factors, boostrapped.
% dates: dates of the discount factors, bootstraped.
% coupon: coupon of the corresponding bond
% cleanPrice: clean price of the bond, no accrual
% startYear: start year of the bond, format datetime
% endYear: end year of the bond, format datetime
% today: today's date
% OUTPUT:
% as_rate: asset swap rate

% We generate a set of dates for six years
t0 = today;
fixedPaymentDates = datetime(startYear:endYear,3,31);

% check if business day
fixedPaymentDates=check_busday(fixedPaymentDates);

% 29-Mar-2024 holy friday
fixedPaymentDates(3)=fixedPaymentDates(3)-1;

% We decided not to compute accrual, as for the scope of this function
% it is irrelevant. In fact, the formula for the asset swap rate considers the
% difference between the IB bond price and the non defaultable one, so the
% accrual terms get canceled out.
ACT_360 = 2;
C_30_360=6;
deltas = [
    yearfrac(t0, fixedPaymentDates(2), C_30_360);
    yearfrac(fixedPaymentDates(2:end-1),fixedPaymentDates(3:end),C_30_360)';
];

% We compute the needed discount factors interpolating on the boostrap
discounts_fixed = interpolation(discounts, dates, fixedPaymentDates(2:end));

% We price the IB bond, with a face value of 100
IB_price = coupon*sum(deltas.*discounts_fixed')+discounts_fixed(end);
IB_price = 100*IB_price;

% We compute the BPV

start_date = datetime(2023,3,31); % Start from 31-Mar-2023
num_dates = 21; % Number of dates to generate
interval = calmonths(3); % 3-month interval

floatPaymentDates = start_date + (0:num_dates-1) * interval;

% check business days
floatPaymentDates=check_busday(floatPaymentDates);

% 29-Mar-2024 holy friday
floatPaymentDates(5)=floatPaymentDates(5)-1;

deltas_f = [
    yearfrac(t0, floatPaymentDates(1), ACT_360);
    yearfrac(floatPaymentDates(1:end-1), floatPaymentDates(2:end), ACT_360)';
];

discounts_float = interpolation(discounts, dates, floatPaymentDates);

BPV = 100*sum(deltas_f.*discounts_float');

% Asset swap rate
as_rate = (IB_price-cleanPrice)/BPV;
end