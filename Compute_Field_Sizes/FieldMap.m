function E = FieldMap(x_E, x_p, a, p, box, xi)

% Compute the field at specific points.
%
% INPUTS
% x_E = (N_E-by-3) field positions
% x_p = (N-by-3) particle positions
% a = (N-by-1) particle radii
% p = (N_p-by-3) particle dipoles
% box = (1-by-3) box dimensions
% xi = (scalar) Ewald parameter
%
% OUTPUTS
% E = (N_E-by-3) field at each x_E

% Config parameters
N_p = size(x_p,1); % number of particles
N_E = size(x_E,1); % number of field points

% Parameters for the spectral Ewald method
tol = 10^-3; % error tolerance
r_c = sqrt(-log(tol))/xi; % real space cutoff radius

% Throw a warning if the real space cutoff is larger than half the box size in any dimension
if any(r_c > box/2)
    error('Real space cutoff (r_c = %.3f) larger than half the box size (%.3f, %.3f, %.3f).',r_c,box(1),box(2),box(3))
end

% Get the unique radii a_uniq
a_uniq = uniquetol(a);

% Build the real space table
r_table = [0.5:tol:r_c+10*tol]'; % separation values for real space table
[Ep_perp, Ep_para] = RealSpaceTable(r_table, a_uniq, xi);

% Compute the neighbor list.
[cell,Ncell] = CellList([x_p;x_E],box,r_c); % cell list for both particle and probe positions
type = [ones(N_p,1);2*ones(N_E,1)];
[n1,n2] = NeighborList_Types(cell,Ncell,type,[2,1]); % neighbor list with only probe/particle pairs
n1 = n1 - N_p; % n1 = 1 now corresponds to x_E(1,:)

E = -ComputeField(x_E, x_p, a, p, box, n1, n2, xi, r_table, Ep_perp, Ep_para, tol);

end