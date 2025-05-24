function [price] = digitalPriceBlk(mkt_strikes, mkt_vols, strike, F0, B, T, flag)

if flag == "BLK"
    sigma = interp1(mkt_strikes, mkt_vols, strike, 'spline');
    d2 = log(F0/strike)/(sigma*sqrt(T))-0.5*sigma*sqrt(T);

    price = B*normcdf(d2);
elseif flag == "Smile"
    impvol = @(K) interp1(mkt_strikes, mkt_vols, K, 'spline');
    
    d1 = @(K) (log(F0 / K) + (0.5 * impvol(K)^2) * T) / (impvol(K) * sqrt(T));
    d2 = @(K) d1(K) - impvol(K) * sqrt(T);
    
    vega = @(K) F0*exp(-d1(K).^2/2)/sqrt(2*pi)*sqrt(T);
    
    % We take the derivative of the smile, approximating it with a finite difference.
    % We need to find the two closest strikes to the one we are interested in.
    i = find(mkt_strikes(1:end-1) <= strike  & mkt_strikes(2:end) >= strike, 1);
    if isempty(i)
        error('Strike not found in the vol surface');
    end
    sigma_diff = (mkt_vols(i+1)-mkt_vols(i))/(mkt_strikes(i+1)-mkt_strikes(i));
    
    price = B*(normcdf(d2(strike))-sigma_diff*vega(strike));
end

end