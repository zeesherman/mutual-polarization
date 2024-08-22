function E_0 = ComputeField(x, a, p, eps_p, box, n1, n2, xi, r_table, Ep_perp, Ep_para, tol)

% Computes the external field at the particle positions consistent with
% their induced dipoles.
%
% INPUTS
% x = (N-by-3) particle positions
% a = (N-by-3) particle radii
% p = (N-by-3) particle dipoles
% eps_p = (N-by-1) particle permittivities
% box = (1-by-3) box dimensions
% n1, n2 = (Nb-by-1) neighbor list
% xi = (scalar) Ewald splitting parameter
% r_table = (n-by-1) separation values in table
% Ep_perp = (n-by-m) I-rr component of the field/dipole coupling; m=M*(M+1)/2 for M unique radii
% Ep_para = (n-by-m) rr component of the field/dipole coupling
% tol = (scalar) error tolerance
%
% OUTPUTS
% E0 = (N-by-3) external field at particle centers

% Parameters of the spectral Ewald method
r_c = sqrt(-log(tol))/xi; % real space cutoff radius
q_c = 2*xi*sqrt(-log(tol)); % wave space cutoff
N_grid = ceil(1+box*q_c/pi); % number of grid nodes in each dimension
h = box./N_grid; % grid spacing in each dimension
%P = ceil(-2*log(tol)/pi); % number of grid nodes over which to support the Gaussians
P = ceil(2/mean(h)-2*log(tol)/pi);
eta = P*(mean(h)*xi).^2/pi; % spectral splitting parameter

% Spread particle dipoles (as charges) to a regular grid
E = Spread(x, a, p, N_grid, h, xi, eta, P);

% Fourier transform the grid
fE = fftshift(fftn(E));

% Wavevectors at each grid point
q_x = (-ceil((N_grid(1)-1)/2):floor((N_grid(1)-1)/2))*2*pi/box(1);
q_y = (-ceil((N_grid(2)-1)/2):floor((N_grid(2)-1)/2))*2*pi/box(2);
q_z = (-ceil((N_grid(3)-1)/2):floor((N_grid(3)-1)/2))*2*pi/box(3);
[Q_X,Q_Y,Q_Z] = ndgrid(q_x,q_y,q_z);
q = cat(4,Q_X,Q_Y,Q_Z);

% Scale the transformed grid
fE_tilde = Scale(fE, q, N_grid, xi, eta);

% Inverse transform the scaled grid
E_tilde = ifftn(ifftshift(fE_tilde));

% Contract the grid to the field at the particles
E_q = Contract(x, a, E_tilde, N_grid, h, xi, eta, P);

% Calculate the real space contribution to the field
E_r = RealSpace(x, a, p, eps_p, box, n1, n2, r_c, r_table, Ep_perp, Ep_para);

% Add the real space and reciprocal space contributions to the field
E_0 = E_q + E_r;

end