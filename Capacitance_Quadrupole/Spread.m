function [HSx,HSy,HSz,HQxx,HQxy,HQxz,HQyy,HQyz] = Spread(x,S,Q,Ngrid,h,xi,eta,P,offset,offsetxyz)

% Spread particle dipoles and quadrupoles to a regular grid.
%
% INPUTS
% x = N-by-3 array of particle positions
% S = N-by-3 array of particle dipoles
% Q = N-by-5 array of particle quadrupoles
% Ngrid = (3 row vector) number of grid nodes in each dimension
% h = (3 row vector) grid spacing in each dimension
% xi = (scalar) Ewald splitting parameter
% eta = (1-by-3) spectral splitting parameter
% P = (scalar) number of grid nodes in each dimension to support Gaussian
%
% OUTPUTS
% HSx = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3) array) x component of the gridded dipoles
% HSy = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3) array) y component of the gridded dipoles
% HSz = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3) array) z component of the gridded dipoles

% Other quantities
N = size(x,1); % number of particles

% Initializations
HSx = zeros(Ngrid);
HSy = zeros(Ngrid);
HSz = zeros(Ngrid);
HQxx = zeros(Ngrid);
HQxy = zeros(Ngrid);
HQxz = zeros(Ngrid);
HQyy = zeros(Ngrid);
HQyz = zeros(Ngrid);

% Express the particle positions in units of grid spacing h
reph = repmat(h,N,1); % copies h into an N-by-3 array
xscaled = x./reph;

% Find the grid node closest to each particle
nearestnode = round(xscaled); % round particle position to nearest node

% Spread the moments of each particle to a grid
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
    HSx(linnode) = HSx(linnode) + Hcoeff*S(n,1);
    HSy(linnode) = HSy(linnode) + Hcoeff*S(n,2);
    HSz(linnode) = HSz(linnode) + Hcoeff*S(n,3);
    HQxx(linnode) = HQxx(linnode) + Hcoeff*Q(n,1);
    HQxy(linnode) = HQxy(linnode) + Hcoeff*Q(n,2);
    HQxz(linnode) = HQxz(linnode) + Hcoeff*Q(n,3);
    HQyy(linnode) = HQyy(linnode) + Hcoeff*Q(n,4);
    HQyz(linnode) = HQyz(linnode) + Hcoeff*Q(n,5);
  
end

end