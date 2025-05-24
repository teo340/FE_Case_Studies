function [dates, swap_rates] = swap_spline(dates, rates)
% We build a complete set of swap rates using a spline interpolation.
% INPUT
% dates: dates of the swaps, format datetime. Read from market data
% rates: rates of the swaps. Read from market data.
% OUTPUT
% dates: complete set of dates, format datetime. Interpolatated
% swap_rates: complete set of swap rates Interpolated

today = datetime('02-Feb-2023'); % settlement date
yf_reduce = yearfrac(today, dates); % year fractions of the market given dates

% The complete set goes from 02-Feb-2024 to 02-Feb-2073.
% Not every date is a business day, so we need to adjust them,
% according to the MODIFIED FOLLOWING convention.
startYear = 2024;
endYear = 2073;
dates = datetime(startYear:endYear,2,2);

% Check if business day
for i = 1:length(dates)
    if weekday(dates(i)) == 1 
        % Sunday - modified following convention
        dates(i) = dates(i) + 1;
    elseif weekday(dates(i)) == 7
        % Saturday - modified following convention
        dates(i) = dates(i) + 2;
    end
end

% Daycount convention: ACT_365, as specificied in the assignment
ACT_365 = 3;
yf = yearfrac(today, dates, ACT_365);
% Third order spline interpolation
swap_rates = spline(yf_reduce, rates, yf);

end     % function swap_spline