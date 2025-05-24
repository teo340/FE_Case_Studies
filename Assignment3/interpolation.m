function [discounts]=interpolation(discounts, dates, target)
% Interpolates the zero rates between the given dates and returns
% the corresponding discount factor. Interpolation is always carried out
% on the zero rates, using a linear interpolator.
%
% INPUT
% discounts: discount factors that contain the target one
% dates: dates in between we interpolate, datenum format
% target: date at which we want to compute the discount factor

% Zero rates are interpolated with yearfrac Convention ACT/365
ACT_365 = 3;

% to compute zero rates, we need the discount factors
% with respect to the settlement date.
today = datenum('02-Feb-2023'); 
zeroRates = -log(discounts)./yearfrac(dates, today, ACT_365);

% Because of how the interp1 MATLAB function works, it's not
% necessary to specify whether to extrapolate or not.
% The function understands what to do itself, based on the values
% of target and the given dates.
% Also, if the target actually coincides with one of the given dates,
% it just returns the corresponding zero rate.

interp_zr = interp1(dates, zeroRates, target, 'linear', zeroRates(end));

% compute discount factor
discounts = exp(-interp_zr.*yearfrac(target, dates(1), ACT_365));

end        % function interpolation