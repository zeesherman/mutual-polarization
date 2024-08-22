function [Ep_perp, Ep_para] = RealSpaceTable(r,a,xi)

% Tabulate the real space contributions as a function of particle separation 
%
% INPUTS
% r = (n-by-1) separation values to tabulate; r must include 0
% a = (1-by-M) distinct radii 
% xi = (scalar) Ewald splitting parameter
%
% OUTPUTS
% Ep_perp = (n-by-m) I-rr component of the field/dipole coupling; m = M*(M+1)/2
% Ep_para = (n-by-m) rr component of the field/dipole coupling

% Display warning if r(1) isn't zero
if r(1) ~= 0
    warning('First entry of real space table is the r->0 limit, but r(1) = %.3g', r(1))
end

% Initialize
n = length(r); % number of separations
M = length(a); % number of radii types
m = M*(M+1)/2; % number of unique radii pairs
Ep_perp = zeros(n,m);
Ep_para = zeros(n,m);
lin = 0; % loop counter

% Loop through m unique radii pairs
for i = 1:M
    for j = i:M
        
        lin = lin + 1;
        a_i = a(i);
        a_j = a(j);
        
        % I-rr component        
        f_1 = 1./(1024*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( 4*r.^5*xi^4 - 4*(a_i+a_j)*r.^4*xi^4 + r.^3*(16*xi^2 + 8*(-4*a_i^2+a_i*a_j-4*a_j^2)*xi^4) ...
            + r.^2*(-12*(a_i+a_j)*xi^2 - 8*(4*a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + 4*a_j^3)*xi^4) ...
            + r*(3 - 12*(a_i^2 - a_i*a_j + a_j^2)*xi^2 - 4*(a_i+a_j)^2*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4) ...
            + 3*(a_i+a_j) + 4*(4*a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + 4*a_j^3)*xi^2 + 4*(a_i+a_j)^3*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4);
        
        f_2 = 1./(1024*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( 4*r.^5*xi^4 + 4*(a_i+a_j)*r.^4*xi^4 + r.^3*(16*xi^2 + 8*(-4*a_i^2+a_i*a_j-4*a_j^2)*xi^4) ...
            + r.^2*(12*(a_i+a_j)*xi^2 + 8*(4*a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + 4*a_j^3)*xi^4) ...
            + r*(3 - 12*(a_i^2 - a_i*a_j + a_j^2)*xi^2 - 4*(a_i+a_j)^2*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4) ...
            - 3*(a_i+a_j) - 4*(4*a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + 4*a_j^3)*xi^2 - 4*(a_i+a_j)^3*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4);
        
        f_3 = 1./(1024*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( -4*r.^5*xi^4 + 4*(-a_i+a_j)*r.^4*xi^4 + r.^3*(-16*xi^2 + 8*(4*a_i^2+a_i*a_j+4*a_j^2)*xi^4) ...
            + r.^2*(12*(-a_i+a_j)*xi^2 - 8*(4*a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - 4*a_j^3)*xi^4) ...
            + r*(-3 + 12*(a_i^2 + a_i*a_j + a_j^2)*xi^2 + 4*(a_i-a_j)^2*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4) ...
            + 3*(a_i-a_j) + 4*(4*a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - 4*a_j^3)*xi^2 + 4*(a_i-a_j)^3*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4);
        
        f_4 = 1./(1024*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( -4*r.^5*xi^4 + 4*(a_i-a_j)*r.^4*xi^4 + r.^3*(-16*xi^2 + 8*(4*a_i^2+a_i*a_j+4*a_j^2)*xi^4) ...
            + r.^2*(12*(a_i-a_j)*xi^2 + 8*(4*a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - 4*a_j^3)*xi^4) ...
            + r*(-3 + 12*(a_i^2 + a_i*a_j + a_j^2)*xi^2 + 4*(a_i-a_j)^2*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4) ...
            + 3*(-a_i+a_j) - 4*(4*a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - 4*a_j^3)*xi^2 - 4*(a_i-a_j)^3*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4);
    
        f_5 = 1./(2048*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( -8*r.^6*xi^6 + 36*r.^4*xi^4*(-1 + 2*(a_i^2+a_j^2)*xi^2) ...
            + 128*(a_i^3+a_j^3)*r.^3*xi^6 + r.^2*(-18*xi^2 + 72*(a_i^2+a_j^2)*xi^4 + 72*(a_i^2-a_j^2)^2*xi^6) ...
            + 3 - 18*(a_i^2+a_j^2)*xi^2 - 36*(a_i^2-a_j^2)^2*xi^4 - 8*(a_i+a_j)^4*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^6);
        
        f_6 = 1./(2048*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( -8*r.^6*xi^6 + 36*r.^4*xi^4*(-1 + 2*(a_i^2+a_j^2)*xi^2) ...
            - 128*(a_i^3+a_j^3)*r.^3*xi^6 + r.^2*(-18*xi^2 + 72*(a_i^2+a_j^2)*xi^4 + 72*(a_i^2-a_j^2)^2*xi^6) ...
            + 3 - 18*(a_i^2+a_j^2)*xi^2 - 36*(a_i^2-a_j^2)^2*xi^4 - 8*(a_i+a_j)^4*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^6);
        
        f_7 = 1./(2048*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( 8*r.^6*xi^6 + 36*r.^4*xi^4*(1 - 2*(a_i^2+a_j^2)*xi^2) ...
            + 128*(a_i^3-a_j^3)*r.^3*xi^6 + r.^2*(18*xi^2 - 72*(a_i^2+a_j^2)*xi^4 - 72*(a_i^2-a_j^2)^2*xi^6) ...
            - 3 + 18*(a_i^2+a_j^2)*xi^2 + 36*(a_i^2-a_j^2)^2*xi^4 + 8*(a_i-a_j)^4*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^6);
        
        f_8 = 1./(2048*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( 8*r.^6*xi^6 + 36*r.^4*xi^4*(1 - 2*(a_i^2+a_j^2)*xi^2) ...
            - 128*(a_i^3-a_j^3)*r.^3*xi^6 + r.^2*(18*xi^2 - 72*(a_i^2+a_j^2)*xi^4 - 72*(a_i^2-a_j^2)^2*xi^6) ...
            - 3 + 18*(a_i^2+a_j^2)*xi^2 + 36*(a_i^2-a_j^2)^2*xi^4 + 8*(a_i-a_j)^4*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^6);
 
        % Regularization for overlapping particles
        f_reg = zeros(n,1);
        
        ind = (r < a_i+a_j) & (r >= a_i-a_j) & (r >= a_j-a_i);
        f_reg(ind) = (a_i + a_j - r(ind)).^3./(128*pi*a_i^3*a_j^3*r(ind).^3) .* ( -r(ind).^3 - 3*(a_i+a_j)*r(ind).^2 ...
            + 3*(a_i^2 - 4*a_i*a_j + a_j^2)*r(ind) + a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + a_j^3 );
        
        ind = (r < a_j-a_i) & (r > a_i-a_j);
        f_reg(ind) = (r(ind).^3 - a_j^3)./(4*pi*a_j^3*r(ind).^3);
        
        ind = (r < a_i-a_j) & (r > a_j-a_i);
        f_reg(ind) = (r(ind).^3-a_i^3)./(4*pi*a_i^3*r(ind).^3);
        
        % Combine
        Ep_perp(:,lin) = f_1.*exp(-(r+a_i+a_j).^2*xi^2) + f_2.*exp(-(r-a_i-a_j).^2*xi^2) + f_3.*exp(-(r-a_i+a_j).^2*xi^2) ...
            + f_4.*exp(-(r+a_i-a_j).^2*xi^2) + f_5.*erfc((r+a_i+a_j)*xi) + f_6.*erfc((r-a_i-a_j)*xi) + f_7.*erfc((r-a_i+a_j)*xi) ...
            + f_8.*erfc((r+a_i-a_j)*xi) + f_reg;
        
        % Handle the r->0 values separately
        Ep_perp(1,lin) = 1/(16*pi^(3/2)*a_i^3*a_j^3*xi^3)*( (1 - 2*a_i^2*xi^2 + 2*a_i*a_j*xi^2 - 2*a_j^2*xi^2)*exp(-(a_i+a_j)^2*xi^2) ...
            + (-1 + 2*a_i^2*xi^2 + 2*a_i*a_j*xi^2 + 2*a_j^2*xi^2)*exp(-(a_i-a_j)^2*xi^2) ) ...
            + 1/(8*pi*a_i^3*a_j^3)*( (a_i^3-a_j^3)*erf((a_i-a_j)*xi) - (a_i^3+a_j^3)*erf((a_i+a_j)*xi) ) + min(a_i^3,a_j^3)/(4*pi*a_i^3*a_j^3);
        
        % rr component
        f_1 = 1./(512*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( 8*r.^5*xi^4 - 8*(a_i+a_j)*r.^4*xi^4 + r.^3*(14*xi^2 - 4*(7*a_i^2 - 4*a_i*a_j + 7*a_j^2)*xi^4) ...
            + r.^2*(-6*(a_i+a_j)*xi^2 - 4*(a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + a_j^3)*xi^4) ...
            + r*(-3 + 12*(a_i^2 - a_i*a_j + a_j^2)*xi^2 + 4*(a_i+a_j)^2*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4) ...
            - 3*(a_i+a_j) - 4*(4*a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + 4*a_j^3)*xi^2 - 4*(a_i+a_j)^3*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4 );
        
        f_2 = 1./(512*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( 8*r.^5*xi^4 + 8*(a_i+a_j)*r.^4*xi^4 + r.^3*(14*xi^2 - 4*(7*a_i^2 - 4*a_i*a_j + 7*a_j^2)*xi^4) ...
            + r.^2*(6*(a_i+a_j)*xi^2 + 4*(a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + a_j^3)*xi^4) ...
            + r*(-3 + 12*(a_i^2 - a_i*a_j + a_j^2)*xi^2 + 4*(a_i+a_j)^2*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4) ...
            + 3*(a_i+a_j) + 4*(4*a_i^3 - 3*a_i^2*a_j - 3*a_i*a_j^2 + 4*a_j^3)*xi^2 + 4*(a_i+a_j)^3*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^4 );
        
        f_3 = 1./(512*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( -8*r.^5*xi^4 + 8*(-a_i+a_j)*r.^4*xi^4 + r.^3*(-14*xi^2 + 4*(7*a_i^2 + 4*a_i*a_j + 7*a_j^2)*xi^4) ...
            + r.^2*(6*(-a_i+a_j)*xi^2 - 4*(a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - a_j^3)*xi^4) ...
            + r*(3 - 12*(a_i^2 + a_i*a_j + a_j^2)*xi^2 - 4*(a_i-a_j)^2*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4) ...
            + 3*(-a_i+a_j) - 4*(4*a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - 4*a_j^3)*xi^2 - 4*(a_i-a_j)^3*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4 );
        
        f_4 = 1./(512*pi^(3/2)*a_i^3*a_j^3*r.^3*xi^5) .* ( -8*r.^5*xi^4 + 8*(a_i-a_j)*r.^4*xi^4 + r.^3*(-14*xi^2 + 4*(7*a_i^2 + 4*a_i*a_j + 7*a_j^2)*xi^4) ...
            + r.^2*(6*(a_i-a_j)*xi^2 + 4*(a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - a_j^3)*xi^4) ...
            + r*(3 - 12*(a_i^2 + a_i*a_j + a_j^2)*xi^2 + 4*(a_i-a_j)^2*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4) ...
            + 3*(a_i-a_j) + 4*(4*a_i^3 + 3*a_i^2*a_j - 3*a_i*a_j^2 - 4*a_j^3)*xi^2 + 4*(a_i-a_j)^3*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^4 );
        
        f_5 = 1./(1024*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( -16*r.^6*xi^6 + 36*r.^4*xi^4*(-1 + 2*(a_i^2+a_j^2)*xi^2) + 64*(a_i^3+a_j^3)*r.^3*xi^6 ...
            - 3 + 18*(a_i^2+a_j^2)*xi^2 + 36*(a_i^2-a_j^2)^2*xi^4 + 8*(a_i+a_j)^4*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^6 );
        
        f_6 = 1./(1024*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( -16*r.^6*xi^6 + 36*r.^4*xi^4*(-1 + 2*(a_i^2+a_j^2)*xi^2) - 64*(a_i^3+a_j^3)*r.^3*xi^6 ...
            - 3 + 18*(a_i^2+a_j^2)*xi^2 + 36*(a_i^2-a_j^2)^2*xi^4 + 8*(a_i+a_j)^4*(a_i^2 - 4*a_i*a_j + a_j^2)*xi^6 );
        
        f_7 = 1./(1024*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( 16*r.^6*xi^6 + 36*r.^4*xi^4*(1 - 2*(a_i^2+a_j^2)*xi^2) + 64*(a_i^3-a_j^3)*r.^3*xi^6 ...
            + 3 - 18*(a_i^2+a_j^2)*xi^2 - 36*(a_i^2-a_j^2)^2*xi^4 - 8*(a_i-a_j)^4*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^6 );
        
        f_8 = 1./(1024*pi*a_i^3*a_j^3*r.^3*xi^6) .* ( 16*r.^6*xi^6 + 36*r.^4*xi^4*(1 - 2*(a_i^2+a_j^2)*xi^2) + 64*(-a_i^3+a_j^3)*r.^3*xi^6 ...
            + 3 - 18*(a_i^2+a_j^2)*xi^2 - 36*(a_i^2-a_j^2)^2*xi^4 - 8*(a_i-a_j)^4*(a_i^2 + 4*a_i*a_j + a_j^2)*xi^6 );
        
        % Regularization for overlapping particles
        f_reg = zeros(n,1);
        
        ind = (r < a_i+a_j) & (r >= a_i-a_j) & (r >= a_j-a_i);
        f_reg(ind) = -1./(64*pi*a_i^3*a_j^3*r(ind).^3) .* ( (a_i+a_j)^4*(a_i^2 - 4*a_i*a_j + a_j^2) - 8*(a_i^3+a_j^3)*r(ind).^3 + 9*(a_i^2+a_j^2)*r(ind).^4 - 2*r(ind).^6 );
        
        ind = (r < a_j-a_i) & (r > a_i-a_j);
        f_reg(ind) = (r(ind).^3+2*a_j^3)./(4*pi*a_j^3*r(ind).^3);
        
        ind = (r < a_i-a_j) & (r > a_j-a_i);
        f_reg(ind) = (r(ind).^3+2*a_i^3)./(4*pi*a_i^3*r(ind).^3);
        
        % Combine
        Ep_para(:,lin) = f_1.*exp(-(r+a_i+a_j).^2*xi^2) + f_2.*exp(-(r-a_i-a_j).^2*xi^2) + f_3.*exp(-(r-a_i+a_j).^2*xi^2) ...
            + f_4.*exp(-(r+a_i-a_j).^2*xi^2) + f_5.*erfc((r+a_i+a_j)*xi) + f_6.*erfc((r-a_i-a_j)*xi) + f_7.*erfc((r-a_i+a_j)*xi) ...
            + f_8.*erfc((r+a_i-a_j)*xi) + f_reg;
        
        % Handle the r->0 values separately. I-rr and rr are the same.
        Ep_para(1,lin) = Ep_perp(1,lin);
        
    end
end

end