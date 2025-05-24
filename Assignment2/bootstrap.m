function [dates, discounts]=bootstrap(datesSet, ratesSet)
% This function bootstraps the discount factor curve
% from the available market data. We consider three different
% kind of instruments: deposits, futures and swaps.
% INPUTS
% datesSet: structure with the dates of the deposits
% ratesSet: structure with the rates of the deposits
% OUTPUT
% dates: dates of the deposits, format datenum
% discounts: discount factors of the deposits
[dates, discounts] = bootstrapDeposits(datesSet, ratesSet);
[dates, discounts] = bootstrapFutures(datesSet, ratesSet, dates, discounts);
[dates, discounts] = bootstrapSwaps(datesSet, ratesSet, dates, discounts);
end     % function bootstrap