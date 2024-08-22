function [S, Q]  = MatrixVectorSolve(x, lambda_p, E0, S_guess, Q_guess, box, p1, p2, Ngrid, h, ...
    P, xi, eta, rc, errortol, offset, offsetxyz, r_vals, field_dip_1, field_dip_2, field_quad_1, ...
    field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4)

% Solves the matrix/vector equation:
% [ M_ES  M_EQ ]   [S]   [ E0 ]
% [ M_GS  M_GQ ] * [Q] = [  0 ]
%
% INPUTS
% x = (N-by-3) particle positions
% lambdap = (N-by-1) particle conductivities
% E0 = (1-by-3) external field
% S_guess = (N-by-3) guess for particle dipoles
% Q_guess = (N-by-5) guess for particle quadrupoles
% box = (1-by-3) dimensions of simulation box
% p1, p2 = (Nnb-by-1) neighbor list
% Ngrid = (1-by3) number of nodes in each dimension to grid simulation box
% h = (1-by-3) grid spacing in each dimension
% P = (scalar) number of grid points in each dimension over which to support Gaussians
% xi = (scalar) Ewald splitting parameter
% eta = (1-by-3) spectral splitting parameter
% rc = (scalar) real space cutoff radius
% errortol = (scalar) error tolerance
% offset = (P^3-by-3) indices of the P^3 surrounding nodes relative to a center node at [0,0,0]
% offsetxyz = (P^3-by-3) coordinates of the P^3 surrounding nodes relative to a center node at [0,0,0]
% r_vals = (Nr-by-1) separation values in tables
% field_dip_1, ..., grad_quad_4 = (Nr-by-1) real space tables
%
% OUTPUTS
% S = (N-by-3) particle dipoles
% Q = (N-by-5) particle quadrupoles

% Initialize
N = size(x,1); % number of particles 
objective = [repmat(E0.',N,1); zeros(5*N,1)]; % right side of equation 
guess = [reshape(S_guess.', 3*N, 1); reshape(Q_guess.', 5*N, 1)]; % reshape initial guess

% Solve the equation iteratively to the specified error tolerance using GMRES.
[sol,~] = gmres(@fun,objective,10,errortol,100,[],[],guess);

% Extract quantities from the solution
S = reshape(sol(1:3*N),3,N).'; % reshape the 3N-by-1 array into a N-by-3 array of particle dipoles
Q = reshape(sol(3*N+1:8*N),5,N).'; % reshape the 5N-by-1 array into a N-by-5 array of particle Quadrupoles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = fun(input)
    
    % Reshape the 8N column vector into an N-by-3 array of particle dipoles
    % and an N-by-5 array of quadrupoles
    S_fun = reshape(input(1:3*N),3,N).';
    Q_fun = reshape(input(3*N+1:8*N),5,N).'; 
    
    % Do the matrix-vector multiply
    [E0_fun, G0_fun] = MatrixVectorMultiply(x, S_fun, Q_fun, lambda_p,box, p1, p2, Ngrid, h, ...
        P, xi, eta, rc, offset, offsetxyz, r_vals, field_dip_1, field_dip_2, field_quad_1, ...
        field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4);
    
    % Reshape the N-by-3 fields and N-by-5 gradients into an 8N-by-1 array
    output = [reshape(E0_fun.',3*N,1); reshape(G0_fun.',5*N,1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end