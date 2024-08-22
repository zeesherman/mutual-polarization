function E = Contract(x, a, E_tilde, N_grid, h, xi, eta, P)

% Contract the gridded field values to particle centers.
%
% INPUTS
% x = (N-by-3) particle positions
% a = (N-by-1) particle radii
% E_tilde = (N_grid(1)-by-N_grid(2)-by-N_grid(3)) gridded fields
% N_grid = (1-by-3) number of grid nodes in each dimension
% h = (1-by-3) grid spacing in each dimension
% xi = (scalar) Ewald splitting parameter
% eta = (scalar) spectral splitting parameter
% P = (1-by-3) Gaussian support (num. grid nodes) in each dimension
%
% OUTPUTS
% E = (N-by-3) wave space contribution to the external field at particle centers

% Number of particles
N = size(x,1);

% Initializations
E = zeros(N,3);

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

% Contract the grid to a field on each particle
for i = 1:N % loop over particles
    
    % Nearest grid node to the current particle
    node_0 = nearestnode(i,:);
    
    % Positions of the P^3 nodes to contract from
    nodes_xyz = node_0.*h + offset_xyz;
    
    % Distances between the spreading nodes and the current particle
    r = nodes_xyz - x(i,:);
    d = vecnorm(r,2,2);
    r_hat = r./d;
    
    % Current particle radius
    a_i = a(i);
    
    % Contraction kernel on the surrounding grid nodes
    c_i = 3./(8*sqrt(2)*pi^(3/2)*a_i^3*sqrt(eta)*xi*d.^2).*( (eta+4*a_i*xi^2*d).*exp(-2*(d+a_i).^2*xi^2/eta) ...
        - (eta-4*a_i*xi^2*d).*exp(-2*(d-a_i).^2*xi^2/eta) ) .* r_hat;

    % Use asymptotic form for small d
    ind = d < 1e-6;
    c_i(ind,:) = 8*sqrt(2)*xi^5*exp(-2*a_i*xi^2/eta)/(pi^(3/2)*eta^(5/2))*r(ind,:);
    
    % Node indices accounting for periodicity
    node = mod(node_0+offset-[1,1,1], N_grid) + [1,1,1];
    
    % Linear index of the nodes
    linnode = node(:,1) + N_grid(1)*(node(:,2)-1) + N_grid(1)*N_grid(2)*(node(:,3)-1);
    
    % Use the trapezoidal rule to sum the surrounding grid nodes.
    E_i = c_i.*E_tilde(linnode); % integrand values
    E(i,:) = prod(h)*sum(E_i,1);

end

end