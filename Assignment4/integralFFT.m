function I = integralFFT(phi, M, dz, queryPoints)
% This function computes the integral of the Lewis formula
% using the Fast Fourier Transform (FFT).
% INPUTS:
% phi: characteristic function of the forward price process
% M: number of points for the FFT
% dz: step size for the FFT
% queryPoints: points at which to evaluate the integral
% OUTPUTS:
% I: vector of integral values at the query points


% FFT parameters
N = 2^M;

z_1 = -(N-1)/2 * dz;
z = z_1:dz:-z_1;

d_xi = 2 * pi / (N * dz);
xi_1 = -(N-1)/2 * d_xi;
xi = xi_1:d_xi:-xi_1;

% integrand function
f = 1 / (2*pi) *  phi(-xi - 1i/2) ./ (xi.^2 + 1/4);
f_tilde = f .* exp(-1i * z_1 * d_xi .* (0:N-1));

FFT = fft(f_tilde);

prefactor = d_xi * exp(-1i * xi_1 * z);
I = prefactor .* FFT;
I = real(I);

% We interpolate the integral in the relevant points
I = interp1(z, I, queryPoints);

end