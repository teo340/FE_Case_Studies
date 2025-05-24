function sigmaShifted = vegaBucket_10y(sigma)

bp=1e-4;                % Define 1 basis point 

% Define transition weights for gradual increase
weights = [0, 1/4, 1/2, 3/4, 1]';

% Initialize output with input values
sigmaShifted = sigma;

% Apply gradual shift to positions 6-10 of the sigma matrix using weights
sigmaShifted(6:10,:) = sigma(6:10,:) + weights*bp;

% Apply full shift to all positions of the sigma matrix from 11 onwards
sigmaShifted(11:end,:) = sigma(11:end,:) + bp;

end