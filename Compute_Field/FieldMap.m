function E = FieldMap(x_E, x_p, p, box, xi)

% Compute the field at specific points.
%
% INPUTS
% x_E = (N_E-by-3) field positions
% x_p = (N_p-by-3) particle positions
% p = (N_p-by-3) particle dipoles
% box = (1-by-3) simulation box
% xi = (scalar) Ewald splitting parameter
%
% OUTPUTS
% E = (N_E-by-3) field at each x_E

% Set numerical parameters
errortol = 10^-3; % error tolerance

% Config parameters
N_p = size(x_p,1); % number of particles
N_E = size(x_E,1); % number of field points

% Build the real space table
r_table = (0.001:0.001:10)'; % separation values for real space table
[~, ~, field_dip_1, field_dip_2, ~, ~] = RealSpaceTable(r_table, xi);
r_table = [0; r_table]; % prepend 0 to r_table

% Parameters for the spectral Ewald method
rc = sqrt(-log(errortol))/xi; % real space cutoff radius
kcut = 2*xi*sqrt(-log(errortol)); % wave space cutoff
Ngrid = ceil(1+box*kcut/pi); % number of grid nodes in each dimension
h = box./Ngrid; % grid spacing in each dimension
P = ceil(-2*log(errortol)/pi); % number of grid nodes over which to support the Gaussians
eta = P*(h*xi).^2/pi; % spectral splitting parameter
[offset, offsetxyz] = PreCalculations(P, h);

% Throw a warning if the real space cutoff is larger than half the box size
% in any dimension
if any(rc > box/2)
    error('Real space cutoff (rc = %.3f) larger than half the box size (%.3f, %.3f, %.3f). Accuracy to specified error tolerance not guaranteed.',rc,box(1),box(2),box(3))
end

% Compute the neighbor list.
[cell,Ncell] = CellList([x_p;x_E],box,rc); % cell list for both particle and probe positions
type = [ones(N_p,1);2*ones(N_E,1)];
[p1,p2] = NeighborList_Types(cell,Ncell,type,[2,1]); % neighbor list with only probe/particle pairs
p1 = p1 - N_p; % p1 = 1 now corresponds to x_E(1,:)

% Compute the field
E = -MagneticField(x_E,x_p,p,box,p1,p2,Ngrid,h,P,xi,eta,rc,offset,offsetxyz,field_dip_1,field_dip_2,r_table);
    
end