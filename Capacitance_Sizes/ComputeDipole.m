function p = ComputeDipole(x, a, eps_p, E0, box, p_guess, n1, n2, xi, r_table, Ep_perp, Ep_para, tol)

% Use GMRES to iteratively solve for the particle dipoles induced by an
% external electric field.
%
% INPUTS
% x = (N-by-3) particle positions
% a = (N-by-1) particle radii
% eps_p = (N-by-1) particle permittivities
% E0 = (1-by-3) external field
% box = (1-by-3) box dimensions
% p_guess = (N-by-3) initial guess for dipoles
% n1, n2 = (Nb-by-1) neighbor list
% xi = (scalar) Ewald splitting parameter
% r_c = (scalar) real space cutoff radius
% r_table = (n-by-1) separation values in table
% Ep_perp = (n-by-m) I-rr component of the field/dipole coupling; m=M*(M+1)/2 for M unique radii
% Ep_para = (n-by-m) rr component of the field/dipole coupling
% tol = (scalar) error tolerance
%
% OUTPUTS
% p = (N-by-3) particle dipoles

% Number of particles 
N = size(x,1);

% Reshape the exteral field and dipoles
E0 = repmat(E0.', N, 1); % 3N-by-1 vector of E0 replicated N times
p_guess = reshape(p_guess.', 3*N, 1); % 3N-by-1 vector of N dipoles

% Solve for the dipoles
restart = min(3*N, 10); % number of inner iterations is 10 unless the system is small
maxit = min(3*N, 100); % number of outer iterations is 100 unless the system is small
[p,~] = gmres(@fun, E0, restart, tol, maxit, [], [], p_guess);

% Rehape dipoles back to N-by-3
p = reshape(p, 3, N).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E0_prime] = fun(p_prime)
    
    % Reshape 3N-by-1 vector into N-by-3 array of particle dipoles
    p_prime = reshape(p_prime,3,N).';
    
    % Compute the external field at the particle center
    E0_prime = ComputeField(x, a, p_prime, eps_p, box, n1, n2, xi, r_table, Ep_perp, Ep_para, tol);
    
    % Reshape N-by-3 array of external fields into 3N-by-1
    E0_prime = reshape(E0_prime.',3*N,1);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end