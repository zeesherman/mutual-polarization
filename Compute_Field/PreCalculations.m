function [offset,offsetxyz] = PreCalculations(P,h)
%function [offset,offsetxyz,A,C,kmax] = PreCalculations(P,h)

% Performs certain calculations once before simulation begins.
%
% INPUTS
% P = number of grid nodes over which to support the Gaussians
% h = grid spacing in each dimension
% xi = Ewald splitting parameter
%
% OUTPUTS
% offset = indices of the P^3 surrounding nodes relative to a center node
%          at [0,0,0]
% offsetxyz = coordinates of the P^3 surrounding nodes relative to a center
%             node at [0,0,0]

% Calculate the indices and coordinates of the P^3 nodes surrounding a node
linnodes = (0:P^3-1)';
offsetx = mod(linnodes,P) - floor((P-1)/2);
offsety = mod(floor(linnodes/P),P) - floor((P-1)/2);
offsetz = floor(linnodes/P^2) - floor((P-1)/2);
offset = [offsetx,offsety,offsetz];
offsetxyz = offset.*repmat(h,P^3,1);

% Load two body coefficients
% A = [];
% C = [];
% kmax = [];
% load('TwoBodyCoefficients_chiInf_kmax200.mat')

end