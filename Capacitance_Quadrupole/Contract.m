function [Ek,Gk] = Contract(x,Ngrid,h,xi,eta,P,HEtilde,HGtilde,offset,offsetxyz)

% Contract the gridded values to the particle centers.
%
% INPUTS
% x = N-by-3 array of particle positions
% Ngrid = row vector of number of grid nodes in each dimension
% h = row vector containing grid spacing in each dimension
% xi = Ewald splitting parameter
% eta = (1-by-3) spectral splitting parameter
% P = number of grid nodes in each dimension to support Gaussian
% HEtilde = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-3 array) gridded field
% fHGtilde = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-5 array) grided gradient
% offset = indices of the P^3 surrounding nodes relative to a center node
%          at [0,0,0]
% offsetxyz = coordinates of the P^3 surrounding nodes relative to a center
%             node at [0,0,0]
% OUTPUTS
% Ek = (N-by-3 array) reciprocal contribution to the field
% Gk = (N-by-5 array) reciprocal contribution to the field gradient

% Other quantities
N = size(x,1); % number of particles

% Initializations
Ek = zeros(N,3);
Gk = zeros(N,5);

% Express the particle positions in units of grid spacing h
reph = repmat(h,N,1); % copies h into an N-by-3 array
xscaled = x./reph;

% Find the grid node closest to each particle
nearestnode = round(xscaled); % round particle position to nearest node

% Reshape the grid array and isolate each component
HEtilde = reshape(HEtilde,[prod(Ngrid),3]);
HExtilde = HEtilde(:,1);
HEytilde = HEtilde(:,2);
HEztilde = HEtilde(:,3);
HGtilde = reshape(HGtilde,[prod(Ngrid),5]);
HGxxtilde = HGtilde(:,1);
HGxytilde = HGtilde(:,2);
HGxztilde = HGtilde(:,3);
HGyytilde = HGtilde(:,4);
HGyztilde = HGtilde(:,5);

% Contract the grid to a field and gradient on each particle
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
    
    % Accumulate each grid node's contribution to particle field and field gradient
    Ek(n,:) = [sum(Hcoeff.*HExtilde(linnode)),sum(Hcoeff.*HEytilde(linnode)),sum(Hcoeff.*HEztilde(linnode))];
    Gk(n,:) = [sum(Hcoeff.*HGxxtilde(linnode)),sum(Hcoeff.*HGxytilde(linnode)),sum(Hcoeff.*HGxztilde(linnode)),...
               sum(Hcoeff.*HGyytilde(linnode)),sum(Hcoeff.*HGyztilde(linnode))];

end

end