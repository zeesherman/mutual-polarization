function [C,S] = Capacitance(x,lambdap,box,Sguess,xi,errortol,dipoletable1,dipoletable2,rtable,offset,offsetxyz)

% Compute the capacitance of a configuration of particles.
%
% INPUTS
% x = (N-by-3) particle positions
% lambdap = (N-by-1) particle conductivities
% box = (1-by-3) box dimensions
% Sguess = (N-by-3) initial guess for particle dipoles
% xi = (scalar) Ewald parameter
% errortol = (scalar) error tolerance
% dipoletable1 = (Ntable vec) I-rr component of the field/dipole real space coupling
% dipoletable2 = (Ntable vec) rr component of the field/dipole real space coupling
% rtable = (Ntable vec) separation values for table entries
% offset = (P^3-by-3) indices of the P^3 surrounding nodes relative to a center node at [0,0,0]
% offsetxyz = (P^3-by-3) coordinates of the P^3 surrounding nodes relative to a center node at [0,0,0]
%
% OUTPUTS
% C = (scalar) effective capacitance
% S = (N-by-3) induced dipoles

% Initializations
H0 = [0,0,1]; % unit magnetic field

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

% Compute the magnetic dipoles.
S  = MagneticDipole(x,lambdap,H0,box,p1,p2,Ngrid,h,P,xi,eta,rc,...
               Sguess,offset,offsetxyz,dipoletable1,dipoletable2,rtable,errortol);

% Compute the capcitance from the average dipole
C = mean(S(:,3));

end