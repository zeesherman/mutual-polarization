function [E, G] = MatrixVectorMultiply(x_E, x_p, S, Q, box, p1, p2, Ngrid, h, P, ...
    xi, eta, rc, offset, offsetxyz, r_vals, field_dip_1, field_dip_2, field_quad_1, ...
    field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4)

% Computes the field and field gradient at probe points due to due to
% dipoles and quadrupoles.
%
% INPUTS
% x_E = (N_E-by-3) field positions
% x_p = (N_E-by-3) particle positions
% S = (N-by-3 array) particle dipoles
% Q = (N-by-5 array) particle quadrupoles
% box = (1-by-3) dimensions of simulation box
% p1, p2 = (Nnb-by-1) neighbor list
% Ngrid = (1-by3) number of nodes in each dimension to grid simulation box
% h = (1-by-3) grid spacing in each dimension
% P = (scalar) number of grid points in each dimension over which to support Gaussians
% xi = (scalar) Ewald splitting parameter
% eta = (1-by-3) spectral splitting parameter
% rc = (scalar) real space cutoff radius
% offset = (P^3-by-3) indices of the P^3 surrounding nodes relative to a center node at [0,0,0]
% offsetxyz = (P^3-by-3) coordinates of the P^3 surrounding nodes relative to a center node at [0,0,0]
% r_vals = (Nr-by-1) separation values in tables
% field_dip_1, ..., grad_quad_4 = (Nr-by-1) real space tables
%
% OUTPUTS
% E = (N_E-by-3 array) field at probe points
% G = (N_E-by-5 array) field gradient at probe points

% Spread particle dipoles and quadrupoles (as vectors and tensors) to a regular grid
[HSx, HSy, HSz, HQxx, HQxy, HQxz, HQyy, HQyz] = Spread(x_p, S, Q, Ngrid, h, xi, eta, P, offset, offsetxyz);

% Perform a Fourier transform on each component of the grid
fHSx = fftshift(fftn(HSx));
fHSy = fftshift(fftn(HSy));
fHSz = fftshift(fftn(HSz));
fHQxx = fftshift(fftn(HQxx));
fHQxy = fftshift(fftn(HQxy));
fHQxz = fftshift(fftn(HQxz));
fHQyy = fftshift(fftn(HQyy));
fHQyz = fftshift(fftn(HQyz));
fHS = cat(4, fHSx, fHSy, fHSz);
fHQ = cat(4, fHQxx, fHQxy, fHQxz, fHQyy, fHQyz);

% Wave vectors corresponding to the reciprocal grid of the Fourier
% transform
kx = (-ceil((Ngrid(1)-1)/2):floor((Ngrid(1)-1)/2))*2*pi/box(1);
ky = (-ceil((Ngrid(2)-1)/2):floor((Ngrid(2)-1)/2))*2*pi/box(2);
kz = (-ceil((Ngrid(3)-1)/2):floor((Ngrid(3)-1)/2))*2*pi/box(3);
[KX,KY,KZ] = ndgrid(kx,ky,kz);
k = cat(4,KX,KY,KZ);

% Scale the transformed grid
[fHEtilde, fHGtilde] = Scale(fHS, fHQ, k, Ngrid, xi, eta);

% Invert each component of the transformed grid
HEtildex = ifftn(ifftshift(fHEtilde(:,:,:,1)));
HEtildey = ifftn(ifftshift(fHEtilde(:,:,:,2)));
HEtildez = ifftn(ifftshift(fHEtilde(:,:,:,3)));
HGtildexx = ifftn(ifftshift(fHGtilde(:,:,:,1)));
HGtildexy = ifftn(ifftshift(fHGtilde(:,:,:,2)));
HGtildexz = ifftn(ifftshift(fHGtilde(:,:,:,3)));
HGtildeyy = ifftn(ifftshift(fHGtilde(:,:,:,4)));
HGtildeyz = ifftn(ifftshift(fHGtilde(:,:,:,5)));
HEtilde = cat(4, HEtildex, HEtildey, HEtildez);
HGtilde = cat(4, HGtildexx, HGtildexy, HGtildexz, HGtildeyy, HGtildeyz);

% Contract the grid to the field at the particles
[Ek,Gk] = Contract(x_E, Ngrid, h, xi, eta, P, HEtilde, HGtilde, offset, offsetxyz);

% Calculate the real space contribution to the field
[Er,Gr] = RealSpace(x_E, x_p, S, Q, box, p1, p2, rc, r_vals, field_dip_1, field_dip_2, ...
    field_quad_1, field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4);

% Add the real space and reciprocal space contributions to the field.
E = Ek + Er;
G = Gk + Gr;

end