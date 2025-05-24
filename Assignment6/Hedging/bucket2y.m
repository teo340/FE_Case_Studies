function ratesShift = bucket2y(ratesSet)
    bp=1e-4;                % Define 1 basis point 
    ratesShift = ratesSet;  % Initialize output with input values
    
    % Apply full shift to deposit rates
    ratesShift.depos = ratesSet.depos + bp;
    
    % Apply full shift to futures rates
    ratesShift.futures = ratesSet.futures + bp;

    weights = [1, 2/3, 1/3, 0]';

    ratesShift.swaps(1:4,:) = ratesSet.swaps(1:4,:) + weights*bp;
end
