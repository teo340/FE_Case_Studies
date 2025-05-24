function ratesShift = bucket_10y(ratesSet)
    bp=1e-4;                % Define 1 basis point 
    ratesShift = ratesSet;  % Initialize output with input values
    weights = [0, 1/4, 1/2, 3/4, 1]'; % Define transition weights
    
    % Apply gradual shift to positions 5-9 using weights
    ratesShift.swaps(5:9,:) = ratesSet.swaps(5:9,:) + weights*bp;
    
    % Apply full shift to all positions from 10 onwards
    ratesShift.swaps(10:end,:) = ratesSet.swaps(10:end,:) + bp;
end
