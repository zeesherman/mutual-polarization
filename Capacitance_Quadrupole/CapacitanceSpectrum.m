function C = CapacitanceSpectrum(x, box, eps_p, xi)

% Compute the capacitance spectrum.
%
% INPUTS
% x = (N-by-3-by-N_frames) particle positions at each frame
% box = (1-by-3) box dimensions
% eps_p = (N-by-N_k) particle dielectric functions at each wavevector
% xi = (scalar) Ewald parameter
%
% OUTPUTS
% C = (1-by-N_k) capacitance at each wavevector

% Set numerical parameters
errortol = 10^-3; % error tolerance

% Config parameters
N = size(x,1); % number of particles
N_frames = size(x,3); % number of frames
N_k = size(eps_p,2); % number of wavevectors
eta = 4*pi*N/(3*prod(box)); % volume fraction

% Build the real space table
r_table = [0:0.0001:10]'; % separation values for real space table
[~, ~, ~, field_dip_1, field_dip_2, ~, ~, field_quad_1, field_quad_2, ...
    field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4] = RealSpaceTable(r_table, xi);

% Pre-calculations
kcut = 2*xi*sqrt(-log(errortol)); % wave space cutoff
Ngrid = ceil(1+box*kcut/pi); % number of grid nodes in each dimension
h = box./Ngrid; % grid spacing in each dimension
P = ceil(-2*log(errortol)/pi); % number of grid nodes over which to support the Gaussians
[offset, offsetxyz] = PreCalculations(P, h);

% Initializations
C = zeros(N_frames, N_k); % capacitance
Q_guess = zeros(N,5); 

% Loop through configurations
for j = 1:N_frames
    
    % Loop through wavevectors
    parfor i = 1:N_k

        % Construct initial guess
        beta = (eps_p(:,i)-1)./(eps_p(:,i)+2);
        p_guess = [zeros(N,2), 4*pi*beta./(1-beta*eta)].';
        
        % Compute the capacitance
        [C(j,i),~,~] = Capacitance(x, eps_p(:,i), box, p_guess, Q_guess, xi, errortol, offset, offsetxyz, r_table, ...
            field_dip_1, field_dip_2, field_quad_1, field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4);
        
    end
end

% Average the capcitance over configurations
C = mean(C,1);

end