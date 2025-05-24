function ratesShift = bucket5y(ratesSet)
    bp=1e-4;                % Define 1 basis point 
    ratesShift = ratesSet;  % Initialize output with input values
    
    % Define increasing transition weights for position 1-4
    weights1 = [0, 1/3, 2/3, 1]';
    
    % Apply gradual increasing shift to positions 1-4
    ratesShift.swaps(1:4,:) = ratesSet.swaps(1:4,:) + weights1*bp;
    
    % Define decreasing transition weights for positions 5-9
    weights2 = [4/5, 3/5, 2/5, 1/5, 0]';
    
    % Apply gradual decreasing shift to positions 5-9
    ratesShift.swaps(5:9,:) = ratesSet.swaps(5:9,:) + weights2*bp;
end
