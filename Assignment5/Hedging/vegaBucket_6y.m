function sigmaShifted = vegaBucket_6y(sigma)
bp=1e-4;                % Define 1 basis point 

% Define transition weights for gradual decrease
weights = [1, 3/4, 1/2, 1/4, 0]';

% Initialize output with input values
sigmaShifted = sigma;

% Apply full shift to positions 1-5 of the sigma matrix
sigmaShifted(1:5,:) = sigma(1:5,:) + bp;

% Apply gradual decreasing shift to positions 6-10 of the sigma matrix using weights
sigmaShifted(6:10,:) = sigma(6:10,:) + weights*bp;
end