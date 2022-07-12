function [m]  = MagneticDipole(x,lambdap,H,box,p1,p2,Ngrid,h,P,xi,eta,rc,mold,offset,offsetxyz,Hperp,Hpara,rvals,errortol)

% Uses the conjugate gradient method to iteratively solve for the particle
% dipoles required to give a specified magnetic field in the system.

N = size(x,1); % number of particles 

H = H.'; % transpose the row vector to a column vector
Hrep = repmat(H,N,1); % replicate H to create a 3N-by-1 column vector of magnetic fields

mold = reshape(mold.',3*N,1);

% Solve for the dipoles
restart = min(3*N, 10); % number of inner iterations is 10 unless the system is small
maxit = min(3*N, 100); % number of outer iterations is 100 unless the system is small
[m,~] = gmres(@fun,Hrep,restart,errortol,maxit,[],[],mold);
%m = cgs(@fun,Hrep,errortol,[],[],[],mold);

m = reshape(m,3,N).'; % reshape the 3N column vector into a N-by-3 array of particle dipoles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = fun(mprime)
    m = reshape(mprime,3,N).'; % reshape the 3N column vector into a N-by-3 array of particle dipoles
    Hprime = MagneticField(x,m,lambdap,box,p1,p2,Ngrid,h,P,xi,eta,rc,offset,offsetxyz,Hperp,Hpara,rvals);
    H = reshape(Hprime.',3*N,1); % reshape the N-by-3 array of particle magnetic fields into a 3N column vector
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end