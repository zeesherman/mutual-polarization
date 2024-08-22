function E = Spread(x, a, p, N_grid, h, xi, eta, P)

% Spread particle dipoles as charges to a regular grid.
%
% INPUTS
% x = (N-by-3) particle positions
% a = (N-by-1) particle radii
% p = (N-by-3) particle dipoles
% N_grid = (1-by-3) number of grid nodes in each dimension
% h = (1-by-3) grid spacing in each dimension
% xi = (scalar) Ewald splitting parameter
% eta = (scalar) spectral splitting parameter
% P = (1-by-3) Gaussian support (num. grid nodes) in each dimension
%
% OUTPUTS
% E = (N_grid(1)-by-N_grid(2)-by-N_grid(3)) gridded dipoles

% Number of particles
N = size(x,1);

% Initializations
E = zeros(N_grid);

% Calculate the indices and coordinates of the P^3 nodes surrounding a node
linnodes = (0:P^3-1)';
offset_x = mod(linnodes,P) - floor((P-1)/2);
offset_y = mod(floor(linnodes/P),P) - floor((P-1)/2);
offset_z = floor(linnodes/P^2) - floor((P-1)/2);
offset = [offset_x,offset_y,offset_z];
offset_xyz = offset.*h;

% Express particle positions in units of grid spacing h
xscaled = x./h;

% Find the grid node closest to each particle
nearestnode = round(xscaled); % round particle position to nearest node

% Spread the dipole of each particle to a grid
for j = 1:N % loop over particles
    
    % Nearest grid node to the current particle
    node_0 = nearestnode(j,:);
    
    % Positions of the P^3 nodes to spread onto
    nodes_xyz = node_0.*h + offset_xyz;
    
    % Distances between the spreading nodes and the current particle
    r = nodes_xyz - x(j,:); 
    d = vecnorm(r,2,2);
    r_hat = r./d;
    
    % Current particle radius
    a_j = a(j);
    
    % Spread dipole to the surrounding grid nodes
    E_j = 3./(8*sqrt(2)*pi^(3/2)*a_j^3*sqrt(eta)*xi*d.^2).*( (eta+4*a_j*xi^2*d).*exp(-2*(d+a_j).^2*xi^2/eta) ...
        - (eta-4*a_j*xi^2*d).*exp(-2*(d-a_j).^2*xi^2/eta) ) .* (r_hat*p(j,:).');
    
    % Use asymptotic form for small d
    ind = d < 1e-6;
    E_j(ind) = 8*sqrt(2)*xi^5*exp(-2*a_j*xi^2/eta)/(pi^(3/2)*eta^(5/2))*(r(ind,:)*p(j,:).');
    
    % Node indices accounting for periodicity
    node = mod(node_0+offset-[1,1,1], N_grid) + [1,1,1];
    
    % Linear index of the nodes
    linnode = node(:,1) + N_grid(1)*(node(:,2)-1) + N_grid(1)*N_grid(2)*(node(:,3)-1);
    
    % Accumulate current particle's contribution to grid values
    E(linnode) = E(linnode) + E_j;
  
end

end