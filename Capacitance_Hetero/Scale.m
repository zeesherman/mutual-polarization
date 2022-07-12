function [fHtilde] = Scale(fH,k,Ngrid,xi,eta)

% Scale the gridded values for the field.
%
% INPUTS
% fH = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3 array) Fourier transform of dipoles on the reciprocal grid
% k = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3 array) wave vector (4th dimension) 
%     corresponding to each reciprocal grid node (first 3 dimensions)
% Ngrid = (3 row vector) number of reciprocal grid nodes in each dimension
% xi = Ewald splitting parameter
% eta = (3 row vector) spectral splitting parameter
%
% OUTPUTS 
% fHtilde = (same size as fH) scaled grid from which the dipoles can be computed

% magnitude squared of the current wave vector
k2 = sum(k.^2,4);

% index of the k = 0 entry
k0x = ceil((Ngrid(1) - 1)/2) + 1;
k0y = ceil((Ngrid(2) - 1)/2) + 1;
k0z = ceil((Ngrid(3) - 1)/2) + 1;

% k unit vectors
khat = k./repmat(sqrt(k2),1,1,1,3);
khat(k0x,k0y,k0z,:) = 0;

% Dot product of eta and k^2 on the grid
etak2 = dot(repmat(reshape((1-eta),1,1,1,3),size(k(:,:,:,1))),k.^2,4);

% scale the Fourier transform of the grid
Htildecoeff = 9*pi./(2*sqrt(k2)).*besselj(1+1/2,sqrt(k2)).^2.*exp(-etak2/(4*xi^2))./k2;
Htildecoeff(k0x,k0y,k0z) = 0;
fHtilde = repmat(Htildecoeff.*dot(khat,fH,4),[1,1,1,3]).*khat;

end