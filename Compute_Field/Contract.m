function [Hk] = Contract(x,Ngrid,h,xi,eta,P,Htilde,offset,offsetxyz)

% Contract the gridded values to the external field at the particle
% centers.
%
% INPUTS
% x = N-by-3 array of particle positions
% Ngrid = row vector of number of grid nodes in each dimension
% h = row vector containing grid spacing in each dimension
% xi = Ewald splitting parameter
% eta = spectral splitting parameter
% P = number of grid nodes in each dimension to support Gaussian
% Htilde = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3 array) real space grid from which the particle field can be calculated
%
% OUTPUTS
% Hk = (N-by-3 array) reciprocal contribution to the magnetic field

% Other quantities
N = size(x,1); % number of particles

% Initializations
Hk = zeros(N,3);

% Express the particle positions in units of grid spacing h
reph = repmat(h,N,1); % copies h into an N-by-3 array
xscaled = x./reph;

% Find the grid node closest to each particle
nearestnode = round(xscaled); % round particle position to nearest node

% Reshape the grid array and isolate each component
Htilde = reshape(Htilde,[prod(Ngrid),3]);
Hxtilde = Htilde(:,1);
Hytilde = Htilde(:,2);
Hztilde = Htilde(:,3);

% Contract the grid to a field on each particle
for n = 1:N % loop over particles
    
    % Nearest grid node to the current particle
    node0 = nearestnode(n,:);
    
    % Positions of the P^3 nodes to contract from
    nodesxyz = repmat(node0.*h,P^3,1) + offsetxyz;
    
    % Distances between the spreading nodes and the current particle
    r = nodesxyz - repmat(x(n,:),P^3,1);
    
    % Coefficient for contracting the grid
    Hcoeff = (2*xi^2/pi)^(3/2)*sqrt(1/prod(eta))*exp(-2*xi^2*r.^2*(1./eta'));
    
    % Trapezoidal rule
    Hcoeff = Hcoeff*prod(h);
    
    % Node indices accounting for periodicity
    node = mod(repmat(node0,P^3,1) + offset - repmat([1,1,1],P^3,1),repmat(Ngrid,P^3,1)) + repmat([1,1,1],P^3,1);
    
    % Linear index of the nodes
    linnode = node(:,1) + Ngrid(1)*(node(:,2)-1) + Ngrid(1)*Ngrid(2)*(node(:,3)-1);
    
    % Accumulate each grid node's contribution to particle field
    Hk(n,:) = [sum(Hcoeff.*Hxtilde(linnode)),sum(Hcoeff.*Hytilde(linnode)),sum(Hcoeff.*Hztilde(linnode))];

end

end