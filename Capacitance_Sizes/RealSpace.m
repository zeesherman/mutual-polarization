function E = RealSpace(x, a, p, eps_p, box, n1, n2, r_c, r_table, Ep_perp, Ep_para)

% INPUTS
% x = (N-by-3) particle positions
% a = (N-by-3) particle radii
% p = (N-by-3) particle dipoles
% eps_p = (N-by-1) particle permittivities
% box = (1-by-3) box dimensions
% n1, n2 = (Nb-by-1) neighbor list
% r_c = (scalar) real space cutoff radius
% r_table = (n-by-1) separation values in table
% Ep_perp = (n-by-m) I-rr component of the field/dipole coupling; m=M*(M+1)/2 for M unique radii
% Ep_para = (n-by-m) rr component of the field/dipole coupling
%
% OUTPUTS
% E = (N-by-3 array) real space contribution to the field

% Get the unique radii a_uniq and radii indices for each particle
[a_uniq, ~, a_ind] = uniquetol(a);
M = length(a_uniq);

% Build an array to help index into the real space tables. For particles i
% and j, col_ind(a_ind(i), a_ind(j)) is the corresponding column in Ep_perp
% and Ep_para
col_ind = zeros(M,M);
lin = 0;
for i = 1:M
    for j = i:M
        lin = lin + 1;
        col_ind(i,j) = lin;
        col_ind(j,i) = lin;
    end
end

% Dielectric contribution
E = 3*p./(4*pi*a.^3.*(eps_p-1));

% Self term
lin = col_ind(sub2ind([M,M],a_ind,a_ind));
E = E + p.*Ep_perp(1,lin).';

% Loop through neighbor list
for i = 1:length(n1) 
    
    % Closest image distance between the particles
    r_vec = x(n1(i),:) - x(n2(i),:);
    r_vec = r_vec - box.*fix(2*r_vec./box);
    r = sqrt(r_vec*r_vec');
    r_hat = r_vec./r;
    
    % Compute the local field if the distance is less than the cutoff distance
    if r < r_c
        
        % Neighbor dipole
        p_j = p(n2(i),:);

        % Get the table column
        a_ind_i = a_ind(n1(i));
        a_ind_j = a_ind(n2(i));
        col_ind_ij = col_ind(a_ind_i, a_ind_j);
        
        % Find the entries in the real space tables between which to interpolate.
        ind = find(r_table >= r, 1);
        r_high = r_table(ind);
        r_low = r_table(ind-1);

        % Interpolate between the values in the table.
        A = Ep_perp(ind-1,col_ind_ij)*(r_high-r)/(r_high-r_low) + Ep_perp(ind,col_ind_ij)*(r-r_low)/(r_high-r_low);
        B = Ep_para(ind-1,col_ind_ij)*(r_high-r)/(r_high-r_low) + Ep_para(ind,col_ind_ij)*(r-r_low)/(r_high-r_low);

        % Accumulate the n2(i)th particle's contribution to the n1(i)th particle's field
        E(n1(i),:) = E(n1(i),:) + A*(p_j - (p_j*r_hat')*r_hat) + B*(p_j*r_hat')*r_hat;

    end 
     
end

end