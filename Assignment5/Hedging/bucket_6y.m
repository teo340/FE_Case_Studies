function ratesShift = bucket_6y(ratesSet)
    bp=1e-4;                % Define 1 basis point 
    ratesShift = ratesSet;  % Initialize output with input values
    
    % Define increasing transition weights for position 1-5
    weights1 = [0, 1/4, 1/2, 3/4, 1]';
    
    % Apply gradual increasing shift to positions 1-5
    ratesShift.swaps(1:5,:) = ratesSet.swaps(1:5,:) + weights1*bp;
    
    % Define decreasing transition weights for positions 6-9
    weights2 = [3/4, 1/2, 1/4, 0]';
    
    % Apply gradual decreasing shift to positions 6-9
    ratesShift.swaps(6:9,:) = ratesSet.swaps(6:9,:) + weights2*bp;
end
