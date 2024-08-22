function [C, S, Q] = CapacitanceTensor(x, lambda_p, box, S_guess, Q_guess, xi, errortol, offset, offsetxyz, r_vals, ...
    field_dip_1, field_dip_2, field_quad_1, field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4)

% Compute the capacitance tensor of a configuration of particles.
%
% INPUTS
% x = (N-by-3) particle positions
% lambda_p = (N-by-1) particle conductivities
% box = (1-by-3) box dimensions
% S_guess = (N-by-3) initial guess for particle dipoles
% Q_guess = (N-by-5) initial guess for particle quadrupoles
% xi = (scalar) Ewald parameter
% errortol = (scalar) error tolerance
% offset = (P^3-by-3) indices of the P^3 surrounding nodes relative to a center node at [0,0,0]
% offsetxyz = (P^3-by-3) coordinates of the P^3 surrounding nodes relative to a center node at [0,0,0]
% r_vals = (Nr-by-1) separation values in tables
% field_dip_1, ..., grad_quad_4 = (Nr-by-1) real space tables
%
% OUTPUTS
% C = (3-by-3) effective capacitance
% S = (N-by-3-by-3) S(:,:,i) are the induced dipoles for field in ith direction
% Q = (N-by-5-by-3) Q(:,:,i) are the induced quadrupoles for field in ith direction

% Initializations
E0 = [1,0,0; 0,1,0; 0,0,1]; % unit fields
C = zeros(3,3); % capacitance tensor
S = zeros(size(x,1),3,3); % particle dipoles
Q = zeros(size(x,1),5,3); % particle quadrupoles

% Parameters for the spectral Ewald method
rc = sqrt(-log(errortol))/xi; % real space cutoff radius
kcut = 2*xi*sqrt(-log(errortol)); % wave space cutoff
Ngrid = ceil(1+box*kcut/pi); % number of grid nodes in each dimension
h = box./Ngrid; % grid spacing in each dimension
P = ceil(-2*log(errortol)/pi); % number of grid nodes over which to support the Gaussians
eta = P*(h*xi).^2/pi; % spectral splitting parameter

% Throw a warning if the real space cutoff is larger than half the box size
% in any dimension
if any(rc > box/2)
    error('Real space cutoff (rc = %.3f) larger than half the box size (%.3f, %.3f, %.3f). Accuracy to specified error tolerance not guaranteed.',rc,box(1),box(2),box(3))
end

% Compute the neighbor list.
[cell,Ncell] = CellList(x,box,rc);
[p1,p2] = NeighborList(cell,Ncell);

% Loop through three orthogonal field directions
for i = 1:3
 
    % Compute the induced moments
    [S_i, Q_i] = MatrixVectorSolve(x, lambda_p, E0(i,:), S_guess, Q_guess, box, p1, p2, Ngrid, h, ...
        P, xi, eta, rc, errortol, offset, offsetxyz, r_vals, field_dip_1, field_dip_2, ...
        field_quad_1, field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4);
    S(:,:,i) = S_i;
    Q(:,:,i) = Q_i;
    
    % Compute the capcitance from the average dipole
    C(:,i) = mean(S_i,1);    
    
    % Update initial guesses from a column permutation
    S_guess = S_i(:,[3,1,2]);
    Q_guess = [-Q_i(:,1)-Q_i(:,4), Q_i(:,5), Q_i(:,2), Q_i(:,1), Q_i(:,3)]; 
    
end

end