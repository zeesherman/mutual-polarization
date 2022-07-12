function [Hmag] = MagneticField(x,m,lambdap,box,p1,p2,Ngrid,h,P,xi,eta,rc,offset,offsetxyz,Hperp,Hpara,rvals)

% Computes the magnetic field at the particle positions due to their induced
% magnetic dipoles.
%
% INPUTS
% x = (N-by-3 array) particle positions
% m = (N-by-3 array) particle dipoles
% lambdap = (N-by-1) particle conductivities
% box = (3 row vector) dimensions of simulation box
% p1, p2 = (column vectors) neighbor list pairs
% Ngrid = (3 row vector) number of nodes in each dimension to grid simulation box
% h = (3 row vector) grid spacing in each dimension
% P = number of grid points in each dimension over which to support Gaussians
% xi = Ewald splitting parameter
% eta = (3 row vector) spectral splitting parameter
% rc = real space cutoff radius
% offset = indices of the P^3 surrounding nodes relative to a center node
%          at [0,0,0]
% offsetxyz = coordinates of the P^3 surrounding nodes relative to a center
%             node at [0,0,0]
% Hperp = (column vector) I-rr component of the real space field for each 
%        separation in rvals
% Hpara =  (column vector) rr component of the real space field for each 
%         separation in rvals
% rvals = (column vector) separation values for which to compute the real
%         space contribution
%
% OUTPUTS
% Hmag = (N-by-3 array) magnetic field at each of the particles

% Spread particle dipoles (as dipoles) to a regular grid
[Hx,Hy,Hz] = Spread(x,m,Ngrid,h,xi,eta,P,offset,offsetxyz);

% Perform a Fourier transform on each component of the grid
fHx = fftshift(fftn(Hx));
fHy = fftshift(fftn(Hy));
fHz = fftshift(fftn(Hz));
fH = cat(4,fHx,fHy,fHz);

% Wave vectors corresponding to the reciprocal grid of the Fourier
% transform
kx = (-ceil((Ngrid(1)-1)/2):floor((Ngrid(1)-1)/2))*2*pi/box(1);
ky = (-ceil((Ngrid(2)-1)/2):floor((Ngrid(2)-1)/2))*2*pi/box(2);
kz = (-ceil((Ngrid(3)-1)/2):floor((Ngrid(3)-1)/2))*2*pi/box(3);
[KX,KY,KZ] = ndgrid(kx,ky,kz);
k = cat(4,KX,KY,KZ);

% Scale the transformed grid
[fHtilde] = Scale(fH,k,Ngrid,xi,eta);

% Invert each component of the transformed grid
Htildex = ifftn(ifftshift(fHtilde(:,:,:,1)));
Htildey = ifftn(ifftshift(fHtilde(:,:,:,2)));
Htildez = ifftn(ifftshift(fHtilde(:,:,:,3)));
Htilde = cat(4,Htildex,Htildey,Htildez);

% Contract the grid to the field at the particles
Hk = Contract(x,Ngrid,h,xi,eta,P,Htilde,offset,offsetxyz);

% Calculate the real space contribution to the field
Hr = RealSpace(x,m,lambdap,box,p1,p2,rc,Hperp,Hpara,rvals);

% Add the real space and reciprocal space contributions to the field.
Hmag = Hk + Hr;

end