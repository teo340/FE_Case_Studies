function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)
% INPUT:
%   setDate            : Settlement date for the coupon bond
%   couponPaymentDates : Vector of dates when coupon payments occur
%   fixedRate          : Fixed coupon rate of the bond
%   dates              : Complete set of dates used in the bootstrap
%   discounts          : Vector of discount factors from bootstrapping
%
% OUTPUT:
%   MacD               : Macaulay Duration of the coupon bond, representing the weighted average time 
%                        (in years) to receive the bond's cash flows

    % Set the reference settlement date
    t0 = setDate;

    % Compute the year fractions from the settlement date to each coupon payment date 
    % using day count convention 6 (30/360 European).
    yf = yearfrac(t0, couponPaymentDates, 6);

    % Use the discount factors from the bootstrap that correspond to the 7-year payment dates
    discounts = [0.968301310232409; discounts(13:18)];

    % Calculate the numerator of the Macaulay Duration formula:
    % - For coupon payments: Multiply the fixed coupon rate by the product of year fractions and discount factors, then sum
    % - For the final principal repayment: Add the product of the year fraction for the last coupon payment 
    %   and its corresponding discount factor 
    numerator = fixedRate*sum(yf.*discounts)+yearfrac(t0,couponPaymentDates(end),6)*discounts(end);

    % Calculate the denominator, representing the present value of all cash flows:
    % - Sum the discounted coupon payments and the discounted principal (nominal value 1)
    denominator = fixedRate*sum(discounts)+discounts(end);

    % Compute the Macaulay Duration as the ratio of the numerator to the denominator
    MacD = numerator/denominator;
end