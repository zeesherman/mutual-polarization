function [C, p] = CapacitanceTensor(x, a, eps_p, box, p_guess, xi, tol, r_table, Ep_perp, Ep_para)

% Compute the capacitance tensor of a configuration of particles.
%
% INPUTS
% x = (N-by-3) particle positions
% a = (N-by-1) particle radii
% eps_p = (N-by-1) particle permittivities
% box = (1-by-3) box dimensions
% p_guess = (N-by-3) initial guess for particle dipoles
% xi = (scalar) Ewald parameter
% tol = (scalar) error tolerance
% r_table = (n-by-1) separation values in table
% Ep_perp = (n-by-m) I-rr component of the field/dipole coupling; m=M*(M+1)/2 for M unique radii
% Ep_para = (n-by-m) rr component of the field/dipole coupling
%
% OUTPUTS
% C = (1-by-9) effective capacitance
% p = (N-by-9) induced dipoles

% Initializations
E0 = [1,0,0; 0,1,0; 0,0,1]; % unit fields
C = zeros(1,9); % capacitance tensor
p = zeros(size(x,1),9); % particle dipoles

% Parameters for the spectral Ewald method
r_c = sqrt(-log(tol))/xi; % real space cutoff radius

% Throw a warning if the real space cutoff is larger than half the box sizemin any dimension
if any(r_c > box/2)
    error('Real space cutoff (r_c = %.3f) larger than half the box size (%.3f, %.3f, %.3f).',r_c,box(1),box(2),box(3))
end

% Compute the neighbor list.
[cell, N_cell] = CellList(x, box, r_c);
[n1, n2] = NeighborList(cell, N_cell);

% Loop through three orthogonal field directions
for i = 1:3
    
    % Compute dipoles
    p_i = ComputeDipole(x, a, eps_p, E0(i,:), box, p_guess, n1, n2, xi, r_table, Ep_perp, Ep_para, tol);
    
    % Store dipoles
    col_ind = (3*(i-1)+1):(3*i);
    p(:,col_ind) = p_i;

    % Compute the capcitance from the average dipole
    C(col_ind) = mean(p_i,1);    
    
    % Update initial guess from a column permutation
    p_guess = p_i(:,[3,1,2]);
    
end

end