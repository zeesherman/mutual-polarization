function fE_tilde = Scale(fE, q, N_grid, xi, eta)

% Scale the gridded values. This converts dipoles to fields in wave space.
%
% INPUTS
% fE = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)) Fourier transform of dipole grid
% q = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3) wavevector (4th dimension) at each grid node
% N_grid = (1-by-3) number of grid nodes in each dimension
% xi = (scalar) Ewald splitting parameter
% eta = (scalar) spectral splitting parameter
%
% OUTPUTS 
% fE_tilde = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)) scaled grid

% Scale the grid
q_mag = sqrt(sum(q.^2,4)); % wavevector magnitudes
fE_tilde = exp(-(1-eta)*q_mag.^2/(4*xi^2))./q_mag.^2.*fE;

% Remove the q=0 value
q_0 = ceil((N_grid-1)/2) + 1; % index of q=0
fE_tilde(q_0(1),q_0(2),q_0(3)) = 0;

end