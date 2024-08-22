function [Er, Gr] = RealSpace(x_E, x_p, S, Q, box, p1, p2, rc, r_vals, field_dip_1, field_dip_2, ...
    field_quad_1, field_quad_2, field_quad_3, grad_quad_1, grad_quad_2, grad_quad_3, grad_quad_4)

% INPUTS
% x_E = (N_E-by-3) field positions
% x_p = (N-by-3) particle positions
% S = (N-by-3) particle dipoles
% Q = (N-by-5) particle quadrupoles
% box = (1-by-3) box dimensions
% p1, p2 = (Nnb-by-1) neighbor list
% rc = (scalar) real space cutoff radius
% r_vals = (Nr-by-1) separation values in tables
% field_dip_1 = (Nr-by-1) field/dipole coupling: I-rr component
% field_dip_2 = (Nr-by-1) field/dipole coupling: rr component
% field_quad_1 = (Nr-by-1) field/quadrupole or field gradient/dipole coupling: rrr*S component
% field_quad_2 = (Nr-by-1) field/quadrupole or field gradient/dipole coupling: Ir*S+Sr+rS component
% field_quad_3 = (Nr-by-1) field/quadrupole or field gradient/dipole coupling: Ir*S component
% grad_quad_1 = (Nr-by-1) field gradient/quadrupole coupling: Irr:Q component
% grad_quad_2 = (Nr-by-1) field gradient/quadrupole coupling: Q component
% grad_quad_3 = (Nr-by-1) field gradient/quadrupole coupling: Irr:Q+2Q*rr+2rr*Q component
% grad_quad_4 = (Nr-by-1) field gradient/quadrupole coupling: rrrr:Q component
%
% OUTPUTS
% Er = (N-by-3) real space contribution to the external field
% Gr = (N-by-5) real space contribution to the external field gradient

% Intialization
Er = zeros(size(x_E,1),3);
Gr = zeros(size(x_E,1),5);

% Loop through the neighbor list
for i = 1:length(p1) 
    
    % Closest image distance between the particles
    r = x_E(p1(i),:) - x_p(p2(i),:);
    r = r - box.*fix(2*r./box);
    d = sqrt(r*r');
    r = r ./ d;

    % Compute external field and field gradient if the distance is less than the cutoff distance
    if d < rc
        
        % Particle dipoles and quadrupoles
        Sj = S(p2(i),:);
        Qj = Q(p2(i),:);
        
        % Dot products
        rrrdotS = [r(1)^2, r(1)*r(2), r(1)*r(3), r(2)^2, r(2)*r(3)]*(Sj*r');
        IrdotS = [1, 0, 0, 1, 0]*(Sj*r');
        Sr = [Sj(1)*r(1), Sj(1)*r(2), Sj(1)*r(3), Sj(2)*r(2), Sj(2)*r(3)];
        rS = [r(1)*Sj(1), r(1)*Sj(2), r(1)*Sj(3), r(2)*Sj(2), r(2)*Sj(3)];
        Qdotr = [Qj(1)*r(1)+Qj(2)*r(2)+Qj(3)*r(3), ...
                 Qj(2)*r(1)+Qj(4)*r(2)+Qj(5)*r(3), ...
                 Qj(3)*r(1)+Qj(5)*r(2)-(Qj(1)+Qj(4))*r(3)];
        Qdotdotrr = [r(1)^2-r(3)^2, 2*r(1)*r(2), 2*r(1)*r(3), r(2)^2-r(3)^2, 2*r(2)*r(3)]*Qj.';
        IrrdotdotQ = [1, 0, 0, 1, 0]*Qdotdotrr;
        Qdotrr = [Qdotr(1)*r(1), Qdotr(1)*r(2), Qdotr(1)*r(3), Qdotr(2)*r(2), Qdotr(2)*r(3)];
        rrdotQ = [r(1)*Qdotr(1), r(1)*Qdotr(2), r(1)*Qdotr(3), r(2)*Qdotr(2), r(2)*Qdotr(3)];
        rrrrdotdotQ = [r(1)^2, r(1)*r(2), r(1)*r(3), r(2)^2, r(2)*r(3)]*Qdotdotrr;
        
        % Find the entries in the real space tables between which to
        % interpolate.
        interpind = find(floor(r_vals/d),1);
        rhigh = r_vals(interpind);
        rlow = r_vals(interpind-1);

        % Interpolate between the values in the table.
        A = field_dip_1(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            field_dip_1(interpind)*(d - rlow)/(rhigh-rlow);
        B = field_dip_2(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            field_dip_2(interpind)*(d - rlow)/(rhigh-rlow);
        C = field_quad_1(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            field_quad_1(interpind)*(d - rlow)/(rhigh-rlow);
        D = field_quad_2(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            field_quad_2(interpind)*(d - rlow)/(rhigh-rlow);
        E = field_quad_3(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            field_quad_3(interpind)*(d - rlow)/(rhigh-rlow);
        F = grad_quad_1(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            grad_quad_1(interpind)*(d - rlow)/(rhigh-rlow);
        G = grad_quad_2(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            grad_quad_2(interpind)*(d - rlow)/(rhigh-rlow);
        H = grad_quad_3(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            grad_quad_3(interpind)*(d - rlow)/(rhigh-rlow);
        I = grad_quad_4(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            grad_quad_4(interpind)*(d - rlow)/(rhigh-rlow);
        
        % Accumulate the p2(i)th particle's contribution to the p1(i)th
        % particle's field and field gradient
        Er(p1(i),:) = Er(p1(i),:) + A*(Sj-(Sj*r')*r) + B*(Sj*r')*r;
        Er(p1(i),:) = Er(p1(i),:) -1/2*( C*Qdotdotrr*r + 2*D*Qdotr ); % quadrupoles get a factor of 1/2 from the moment expansion
        Gr(p1(i),:) = Gr(p1(i),:) + C*rrrdotS + D*(Sr+rS) + (D+E)*IrdotS;
        Gr(p1(i),:) = Gr(p1(i),:) +1/2*( F*IrrdotdotQ + G*Qj + H*(IrrdotdotQ+2*Qdotrr+2*rrdotQ) + I*rrrrdotdotQ );

    end 
     
end

end