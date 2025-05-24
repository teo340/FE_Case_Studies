function ratesShift = bucket10y(ratesSet)
    bp=1e-4;                % Define 1 basis point 
    ratesShift = ratesSet;  % Initialize output with input values
    weights = [0, 1/5, 2/5, 3/5, 4/5, 1]'; % Define transition weights
    
    % Apply gradual shift to positions 4-9 using weights
    ratesShift.swaps(4:9,:) = ratesSet.swaps(4:9,:) + weights*bp;
    
    % Apply full shift to all positions from 10 onwards
    ratesShift.swaps(10:end,:) = ratesSet.swaps(10:end,:) + bp;
end
