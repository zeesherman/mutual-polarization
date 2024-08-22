function [C, p] = CapacitanceTensorSpectrum(x, box, a, eps_p, xi)

% Compute the capacitance tensor spectrum.
%
% INPUTS
% x = (N-by-3-by-N_frames) particle positions at each frame
% box = (1-by-3) box dimensions
% a = (N-by-1) particle radii
% eps_p = (N-by-N_k) particle dielectric functions at each wavevector
% xi = (scalar) Ewald parameter
%
% OUTPUTS
% C = (N_k-by-9) capacitance tensor at each wavevector; lin components [xx, yx, zx, xy, yy, zy, xz, yz, zz]
% p = (N-by-9-by-N_k-by-N_frames) particle dipoles at each wavevector and frame

% Set numerical parameters
tol = 10^-3; % error tolerance

% Config parameters
N = size(x,1); % number of particles
N_frames = size(x,3); % number of frames
N_k = size(eps_p,2); % number of wavevectors
phi = 4*pi*sum(a.^3)/(3*prod(box)); % volume fraction

% Get the unique radii a_uniq
a_uniq = uniquetol(a);

% Build the real space table
r_table = [0, 1:tol:10]'; % separation values for real space table
[Ep_perp, Ep_para] = RealSpaceTable(r_table, a_uniq, xi);

% Initializations. Wavenumbers need to get the last index because they are
% in the parallel loop
C = zeros(9, N_frames, N_k); % capacitance tensor
p = zeros(N, 9, N_frames, N_k); % particle dipoles

% Loop through configurations
for j = 1:N_frames
    
    % Current configuration
    x_j = x(:,:,j);
    
    % Loop through wavevectors
    parfor i = 1:N_k

        % Construct initial guess
        alpha = (eps_p(:,i)-1)./(eps_p(:,i)+2);
        alpha(isinf(eps_p(:,i))) = 1;
        p_guess = [4*pi*a.^3.*alpha./(1-alpha*phi),zeros(N,2)];
        
        % Compute the capacitance
        [C(:,j,i),p(:,:,j,i)] = CapacitanceTensor(x_j, a, eps_p(:,i), box, p_guess, xi, tol, r_table, Ep_perp, Ep_para);
        
    end
end

% Average the capcitance over configurations
C = squeeze(mean(C,2)).';

% Switch the 3rd and 4th indices; now frames are the last index
p = permute(p,[1,2,4,3]);

end