function [p_ave, p] = DipoleSpectrum(x, box, a, eps_p, E_0, xi, varargin)

% Compute the dipole spectrum.
%
% INPUTS
% x = (N-by-3-by-N_frames) particle positions at each frame
% box = (1-by-3) box dimensions
% a = (N-by-1) particle radii
% eps_p = (N-by-N_k) particle dielectric functions at each wavevector
% E_0 = (1-by-3) incident field polarization
% xi = (scalar) Ewald parameter
% p_guess = (N-by-3-by-N_k-by-N_frames) optional; initial guess for particle dipoles
%
% OUTPUTS
% p_ave = (N_k-by-3) average dipole at each wavevector
% p = (N-by-3-by-N_k-by-N_frames) particle dipoles at each wavevector and frame

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
p_ave = zeros(3, N_frames, N_k); % capacitance tensor
p = zeros(N, 3, N_frames, N_k); % particle dipoles

% Loop through configurations
for j = 1:N_frames
    
    % Current configuration
    x_j = x(:,:,j);
    
    % Loop through wavevectors
    parfor i = 1:N_k
        
        % Construct initial guess if not specified
        if isempty(varargin)
            alpha = (eps_p(:,i)-1)./(eps_p(:,i)+2);
            alpha(isinf(eps_p(:,i))) = 1;
            p_guess = 4*pi*a.^3.*alpha./(1-alpha*phi).*E_0;
        else
            p_guess = varargin{1}(:,:,i,j);
        end
        
        % Compute the capacitance
        [p_ave(:,j,i),p(:,:,j,i)] = Dipole(x_j, a, eps_p(:,i), box, E_0, p_guess, xi, tol, r_table, Ep_perp, Ep_para);
        
    end
end

% Average the dipoles over configurations
p_ave = squeeze(mean(p_ave,2)).';

% Switch the 3rd and 4th indices; now frames are the last index
p = permute(p,[1,2,4,3]);

end