function ratesShift = bucket_2y(ratesSet)
    bp=1e-4;                % Define 1 basis point 
    ratesShift = ratesSet;  % Initialize output with input values
    
    % Apply full shift to deposit rates
    ratesShift.depos = ratesSet.depos + bp;
    
    % Apply full shift to futures rates
    ratesShift.futures = ratesSet.futures + bp;
    
    % Define decreasing transition weights 
    weights = [1, 3/4, 1/2, 1/4, 0]';
    
    % Apply weighted shift to first 5 positions of swap rates
    ratesShift.swaps(1:5,:) = ratesSet.swaps(1:5,:) + weights*bp;
end
