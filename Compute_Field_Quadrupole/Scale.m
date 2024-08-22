function [fHEtilde, fHGtilde] = Scale(fHS, fHQ, k, Ngrid, xi, eta)

% Scale the gridded values in wave space in the matrix/vector multiply.
%
% INPUTS
% fHS = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3 array) Fourier transform of dipoles on the reciprocal grid
% fHQ = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-5 array) Fourier transform of quadrupoles on the reciprocal grid
% k = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3 array) wave vector (4th dimension) 
%     corresponding to each reciprocal grid node (first 3 dimensions)
% Ngrid = (3 row vector) number of reciprocal grid nodes in each dimension
% xi = (scalar) Ewald splitting parameter
% eta = (1-by-3) spectral splitting parameter
%
% OUTPUTS
% fHEtilde = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3 array) scaled field grid
% fHGtilde = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-5 array) scaled field gradient grid

% magnitude squared and magnitude of the current wave vector
k2 = sum(k.^2,4);
kmag = sqrt(k2);

% index of the k = 0 entry
k0x = ceil((Ngrid(1) - 1)/2) + 1;
k0y = ceil((Ngrid(2) - 1)/2) + 1;
k0z = ceil((Ngrid(3) - 1)/2) + 1;

% k unit vectors
khat = k./repmat(kmag,1,1,1,3);
khat(k0x,k0y,k0z,:) = 0;

% Vector versions of the kk-I/3 tensor in the quadrupole shape factor
Qfactor = cat(4, khat(:,:,:,1).*khat(:,:,:,1)-1/3, khat(:,:,:,1).*khat(:,:,:,2), ...
                 khat(:,:,:,1).*khat(:,:,:,3), khat(:,:,:,2).*khat(:,:,:,2)-1/3, ...
                 khat(:,:,:,2).*khat(:,:,:,3)); 
% Qfactor_dot = cat(4, 2*khat(:,:,:,1).*khat(:,:,:,1)+khat(:,:,:,2).*khat(:,:,:,2)-1, ...
%                      2*khat(:,:,:,1).*khat(:,:,:,2), 2*khat(:,:,:,1).*khat(:,:,:,3), ...
%                      2*khat(:,:,:,2).*khat(:,:,:,2)+khat(:,:,:,1).*khat(:,:,:,1)-1, ...
%                      2*khat(:,:,:,2).*khat(:,:,:,3)); % for dotting with Q
                 
Qfactor_dot = cat(4, khat(:,:,:,1).*khat(:,:,:,1)-khat(:,:,:,3).*khat(:,:,:,3), ...
                     2*khat(:,:,:,1).*khat(:,:,:,2), 2*khat(:,:,:,1).*khat(:,:,:,3), ...
                     khat(:,:,:,2).*khat(:,:,:,2)-khat(:,:,:,3).*khat(:,:,:,3), ...
                     2*khat(:,:,:,2).*khat(:,:,:,3)); % for dotting with Q                 

% Dot product of eta and k^2 on the grid
etak2 = dot(repmat(reshape((1-eta),1,1,1,3),size(k(:,:,:,1))),k.^2,4);

% Spherical bessel functions
j1 = sqrt(pi./(2*kmag)).*besselj(1+1/2,kmag);
j2 = sqrt(pi./(2*kmag)).*besselj(2+1/2,kmag);

% Coefficients of the scalings
HEStildecoeff = 9*j1.^2.*exp(-etak2/(4*xi^2))./k2;
HEQtildecoeff = -45/2*1i*j1.*j2.*exp(-etak2/(4*xi^2))./k2;
HGStildecoeff = 45*1i*j2.*j1.*exp(-etak2/(4*xi^2))./k2;
HGQtildecoeff = 225/2*j2.^2.*exp(-etak2/(4*xi^2))./k2;

% Set the coeffcient of the k = 0 term to zero
HEStildecoeff(k0x,k0y,k0z) = 0;
HEQtildecoeff(k0x,k0y,k0z) = 0;
HGStildecoeff(k0x,k0y,k0z) = 0;
HGQtildecoeff(k0x,k0y,k0z) = 0;

% Scale the Fourier transform of the grid
fHEStilde = repmat(HEStildecoeff.*dot(khat,fHS,4),[1,1,1,3]).*khat;
fHEQtilde = repmat(HEQtildecoeff.*dot(Qfactor_dot,fHQ,4),[1,1,1,3]).*khat;
fHGStilde = repmat(HGStildecoeff.*dot(khat,fHS,4),[1,1,1,5]).*Qfactor;
fHGQtilde = repmat(HGQtildecoeff.*dot(Qfactor_dot,fHQ,4),[1,1,1,5]).*Qfactor;

% Combine contributions
fHEtilde = fHEStilde + fHEQtilde;
fHGtilde = fHGStilde + fHGQtilde;

end