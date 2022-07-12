function [p1,p2] = NeighborList(cell,Ncell)

% Compute a list of particle pairs in neighboring cells.
% 
% INPUTS
% cell: index of the cell to which each particle belongs
% Ncell: number of cells in each dimension
%
% OUTPUTS
%
% p1: ID number of the first neighbor pair particle
% p2: ID number of the second neighbor pair particle

% Initializations.
p1 = [];
p2 = [];

% Other quantities.
N = size(cell,1); % number of particles
id = (1:N)'; % particle ID number

% Convert the matrix cell index to a linear index. 
index = cell(:,1) + Ncell(1)*(cell(:,2)-1) + Ncell(1)*Ncell(2)*(cell(:,3)-1);

% Sort the particles according to their cell index.
A = sortrows([index,id,cell]);
index = A(:,1);
id = A(:,2);
cell = A(:,3:5);

% Count the number of particles in each cell.
counts = histc(index,1:prod(Ncell));

% First particle number in each cell.  For example, if cells 1, 2, and 3,
% each contain 3, 4, and 5 particles respectively, particle 1 is the first
% particle in cell 1, particle 4 is the first particle in cell 2, and
% particle 8 is the first particle in cell 3.
cellstart = cumsum(counts);
cellstart = cellstart - counts + 1;

% Loop over all cells.
for i = 1:prod(Ncell)
    
    % Only consider cells that contain particles.
    if counts(i) ~= 0 
        
        % ID number of particles in the current cell.
        pa = id(cellstart(i) - 1 + (1:counts(i)));

        % Search the 27 neighboring (including the current) cells.
        for m = -1:1
            for n = -1:1
                for o = -1:1

                    % Linear cell index of a neighbor cell, taking into
                    % account periodic boundaries.
                    newcell = mod(cell(cellstart(i),:)-1+[m,n,o],Ncell) + 1;
                    newcell = newcell(1) + Ncell(1)*(newcell(2)-1) + Ncell(1)*Ncell(2)*(newcell(3)-1);

                    % Only considering neighboring cells that contain
                    % particles.
                    if counts(newcell) ~= 0
                        % ID number of the particles in the neighbor cell
                        pb = id(cellstart(newcell)-1+(1:counts(newcell)));

                        % Match each particle in the current cell with each
                        % particle in the neighbor cell.
                        [px,py] = meshgrid(pa,pb);
                        px = reshape(px,size(px,1)*size(px,2),1);
                        py = reshape(py,size(py,1)*size(py,2),1);

                        % Update the list of pairs.
                        p1 = [p1;px];
                        p2 = [p2;py];
                    end

                end
            end
        end

    end

end

% Remove pairs of a particle neighboring itself. 
notself = p1 ~= p2;
p1 = p1(notself);
p2 = p2(notself);

end
