function [Hx,Hy,Hz] = Spread(x,m,Ngrid,h,xi,eta,P,offset,offsetxyz)

% Spread particle dipoles to a regular grid.
%
% INPUTS
% x = N-by-3 array of particle positions
% m = N-by-3 array of particle magnetic dipole
% Ngrid = row vector of number of grid nodes in each dimension
% h = row vector containing grid spacing in each dimension
% xi = Ewald splitting parameter
% eta = (3 row vector) spectral splitting parameter
% P = number of grid nodes in each dimension to support Gaussian
%
% OUTPUTS
% Hx = x component of the magnetic dipoles spread to a grid
% Hy = y component of the magnetic dipoles spread to a grid
% Hy = z component of the magnetic dipoles spread to a grid

% Other quantities
N = size(x,1); % number of particles

% Initializations
Hx = zeros(Ngrid);
Hy = zeros(Ngrid);
Hz = zeros(Ngrid);

% Express the particle positions in units of grid spacing h
reph = repmat(h,N,1); % copies h into an N-by-3 array
xscaled = x./reph;

% Find the grid node closest to each particle
nearestnode = round(xscaled); % round particle position to nearest node

% Spread the magnetic dipole of each particle to a grid
for n = 1:N % loop over particles
    
    % Nearest grid node to the current particle
    node0 = nearestnode(n,:);
    
    % Positions of the P^3 nodes to spread onto
    nodesxyz = repmat(node0.*h,P^3,1) + offsetxyz;
    
    % Distances between the spreading nodes and the current particle
    r = nodesxyz - repmat(x(n,:),P^3,1); 

    % Coefficient for spreading the dipoles
    Hcoeff = (2*xi^2/pi)^(3/2)*sqrt(1/prod(eta))*exp(-2*xi^2*r.^2*(1./eta'));

    % Node indices accounting for periodicity
    node = mod(repmat(node0,P^3,1) + offset - repmat([1,1,1],P^3,1),repmat(Ngrid,P^3,1)) + repmat([1,1,1],P^3,1);
    
    % Linear index of the nodes
    linnode = node(:,1) + Ngrid(1)*(node(:,2)-1) + Ngrid(1)*Ngrid(2)*(node(:,3)-1);
    
    % Accumulate current particle's contribution to the spread dipoles
    Hx(linnode) = Hx(linnode) + Hcoeff*m(n,1);
    Hy(linnode) = Hy(linnode) + Hcoeff*m(n,2);
    Hz(linnode) = Hz(linnode) + Hcoeff*m(n,3);
  
end

end