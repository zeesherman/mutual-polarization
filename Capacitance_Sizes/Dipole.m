function [p_ave, p] = Dipole(x, a, eps_p, box, E_0, p_guess, xi, tol, r_table, Ep_perp, Ep_para)

% Compute the mutual dipoles of a configuration of particles.
%
% INPUTS
% x = (N-by-3) particle positions
% a = (N-by-1) particle radii
% eps_p = (N-by-1) particle permittivities
% box = (1-by-3) box dimensions
% E_0 = (1-by-3) incident field polarization
% p_guess = (N-by-3) initial guess for particle dipoles
% xi = (scalar) Ewald parameter
% tol = (scalar) error tolerance
% r_table = (n-by-1) separation values in table
% Ep_perp = (n-by-m) I-rr component of the field/dipole coupling; m=M*(M+1)/2 for M unique radii
% Ep_para = (n-by-m) rr component of the field/dipole coupling
%
% OUTPUTS
% p_ave = (1-by-3) average dipole
% p = (N-by-3) induced dipoles

% Parameters for the spectral Ewald method
r_c = sqrt(-log(tol))/xi; % real space cutoff radius

% Throw a warning if the real space cutoff is larger than half the box sizemin any dimension
if any(r_c > box/2)
    error('Real space cutoff (r_c = %.3f) larger than half the box size (%.3f, %.3f, %.3f).',r_c,box(1),box(2),box(3))
end

% Compute the neighbor list.
[cell, N_cell] = CellList(x, box, r_c);
[n1, n2] = NeighborList(cell, N_cell);
    
% Compute dipoles
p = ComputeDipole(x, a, eps_p, E_0, box, p_guess, n1, n2, xi, r_table, Ep_perp, Ep_para, tol);

% Compute the capcitance from the average dipole
p_ave = mean(p);    

end