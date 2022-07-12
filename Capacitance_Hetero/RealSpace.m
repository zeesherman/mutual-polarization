function [Hr] = RealSpace(x,m,lambdap,box,p1,p2,rc,perp,para,rvals)

% INPUTS
% x = (N-by-3 array) particle positions
% m = (N-by-3 array) particle magnetic dipoles
% lambdap = (N-by-1) particle conductivities
% box = (3 row vector) box dimensions
% p1 = (column vector) first particles in the neighbor pair list
% p2 = (column vector) second particles in the neighbor pair list
% rc = real space cutoff radius
% field_dip_1 = (column vector) I-rr component of the real space field
%                 for each separation in rvals
% field_dip_2 =  (column vector) rr component of the real space field for
%                  each separation in rvals
% rvals = (column vector) separation values for which to compute the real
%         space contribution
%
% OUTPUTS
% Hr = (N-by-3 array) real space contribution to the magnetic field

% Dielectric contribution
Hr = -3./(4*pi*(1-lambdap)).*m;

% Self term
Hr = Hr + m*perp(1);

for i = 1:length(p1) 
    
    % Closest image distance between the particles
    r = x(p1(i),:) - x(p2(i),:);
    r = r - box.*fix(2*r./box);
    d = sqrt(r*r');
    r = r ./ d;

    % Compute potential and force if the distance is less than the cutoff distance
    if d < rc
        
        % Particle dipoles
        mj = m(p2(i),:);
        
        % Find the entries in the real space tables between which to
        % interpolate.
        interpind = find(floor(rvals/d),1);
        rhigh = rvals(interpind);
        rlow = rvals(interpind-1);

        % Interpolate between the values in the table.
        A = perp(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            perp(interpind)*(d - rlow)/(rhigh-rlow);
        B = para(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            para(interpind)*(d - rlow)/(rhigh-rlow);
        
        % Accumulate the p2(i)th particle's contribution to the p1(i)th
        % particle's potential
        Hr(p1(i),:) = Hr(p1(i),:) + A*(mj - (mj*r')*r) + B*(mj*r')*r;

    end 
     
end

end