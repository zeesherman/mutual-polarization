function [cell,Ncell] = CellList(x,box,rmax)

% Bin the particles.
%
% INPUTS
% x: particle positions [ x1 y1 z1; x2 y2 z2; ... ]
% box: a vector containing the x, y, z dimensions of the cell
% rmax: maximum distance over which to look for a pair
%
% OUTPUTS
% cell: index of the cell to which each particle belongs
% Ncell: number of cells in each dimension

% Number of cells in each dimension must be at least 3.  
x = mod(x,box);
Ncell = floor(box/rmax);  
Ncell(Ncell < 3) = 3;

% Replicate the vectors to matrices. 
box = repmat(box,size(x,1),1);
Ncell = repmat(Ncell,size(x,1),1);

% Assign each particle to a cell.
cell = floor(Ncell.*x./box) + 1;

% Recover the vector form of the number of cells.
Ncell = Ncell(1,:);

end
