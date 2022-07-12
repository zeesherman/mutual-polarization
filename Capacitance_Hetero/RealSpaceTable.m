function [pot_charge,pot_dip,field_dip_1,field_dip_2,field_dip_force_1,...
         field_dip_force_2] = RealSpaceTable(r,xi)

% Tabulate the real space contributions as a function of particle separation 
%
% INPUTS
% r = (n vector) separation values to tabulate
% xi = (scalar) Ewald splitting parameter
%
% OUTPUTS
% pot_charge = (n vector) real space contribution to the potential/
%              charge coupling
% pot_dip = (n vector) real space contribution to the potential/
%           dipole or field/charge coupling
% field_dip_1 = (n vector) real space contribution to the I-rr 
%               component of the field/dipole coupling
% field_dip_2 = (n vector) real space contribution to the rr component
%               of the field/dipole coupling
% field_dip_force_1 = (n vector) real space contribution to the
%                     -(mi*mj)r and -((mj*r)mi + (mi*r)mj - 2(mi*r)(mj*r)r) 
%                     component of the field/dipole force
% field_dip_force_2 = (n vector) real space contribution to the
%                     (mi*r)(mj*r)r component of the field/dipole force

%% Potential/Charge Coupling

% Polynomials multiplying the exponetials
exppolyp = -(r+2)./(32*pi^(3/2)*xi*r);
exppolym = -(r-2)./(32*pi^(3/2)*xi*r);
exppoly0 = 1./(16*pi^(3/2)*xi);

% Polynomials multiplying the error functions
erfpolyp = (2*xi^2*(r+2).^2 + 1)./(64*pi*xi^2*r);
erfpolym = (2*xi^2*(r-2).^2 + 1)./(64*pi*xi^2*r);
erfpoly0 = -(2*xi^2*r.^2 + 1)./(32*pi*xi^2*r);

% Regularization for overlapping particles
regpoly = -1./(4*pi*r) + (4-r)./(16*pi);

% Combine the polynomial coefficients, exponentials, and error functions
pot_charge = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi) + (r < 2).*regpoly;

%% Potential/Dipole or Field/Charge coupling

% Polynomials multiplying the exponetials
exppolyp = 1./(256*pi^(3/2)*xi^3*r.^2).*(-6*xi^2*r.^3 - 4*xi^2*r.^2 + (-3+8*xi^2)*r + 2*(1-8*xi^2));
exppolym = 1./(256*pi^(3/2)*xi^3*r.^2).*(-6*xi^2*r.^3 + 4*xi^2*r.^2 + (-3+8*xi^2)*r - 2*(1-8*xi^2));
exppoly0 = 3*(2*r.^2*xi^2+1)./(128*pi^(3/2)*xi^3*r);

% Polynomials multiplying the error functions
erfpolyp = 1./(512*pi*xi^4*r.^2).*(12*xi^4*r.^4 + 32*xi^4*r.^3 + 12*xi^2*r.^2 - 3+64*xi^4);
erfpolym = 1./(512*pi*xi^4*r.^2).*(12*xi^4*r.^4 - 32*xi^4*r.^3 + 12*xi^2*r.^2 - 3+64*xi^4);
erfpoly0 = -3*(4*xi^4*r.^4 + 4*xi^2*r.^2 - 1)./(256*pi*xi^4*r.^2);

% Regularization for overlapping particles
regpoly = -1./(4*pi*r.^2) + r/(8*pi).*(1-3/8*r);

% Combine the polynomial coefficients, exponentials, and error functions
pot_dip = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi) + (r < 2).*regpoly;

%% Field/Dipole coupling: I-rr component

% Polynomials multiplying the exponentials
exppolyp = 1./(1024*pi^(3/2)*xi^5*r.^3).*(4*xi^4*r.^5 - 8*xi^4*r.^4 + 8*xi^2*(2-7*xi^2)*r.^3 - ...
    8*xi^2*(3+2*xi^2)*r.^2 + (3-12*xi^2+32*xi^4)*r + 2*(3+4*xi^2-32*xi^4));
exppolym = 1./(1024*pi^(3/2)*xi^5*r.^3).*(4*xi^4*r.^5 + 8*xi^4*r.^4 + 8*xi^2*(2-7*xi^2)*r.^3 + ...
    8*xi^2*(3+2*xi^2)*r.^2 + (3-12*xi^2+32*xi^4)*r - 2*(3+4*xi^2-32*xi^4));
exppoly0 = 1./(512*pi^(3/2)*xi^5*r.^2).*(-4*xi^4*r.^4 - 8*xi^2*(2-9*xi^2)*r.^2 - 3+36*xi^2);

% Polynomials multiplying the error functions
erfpolyp = 1./(2048*pi*xi^6*r.^3).*(-8*xi^6*r.^6 - 36*xi^4*(1-4*xi^2)*r.^4 + 256*xi^6*r.^3 - ...
    18*xi^2*(1-8*xi^2)*r.^2 + 3-36*xi^2+256*xi^6);
erfpolym = 1./(2048*pi*xi^6*r.^3).*(-8*xi^6*r.^6 - 36*xi^4*(1-4*xi^2)*r.^4 - 256*xi^6*r.^3 - ...
    18*xi^2*(1-8*xi^2)*r.^2 + 3-36*xi^2+256*xi^6);
erfpoly0 = 1./(1024*pi*xi^6*r.^3).*(8*xi^6*r.^6 + 36*xi^4*(1-4*xi^2)*r.^4 + 18*xi^2*(1-8*xi^2)*r.^2 - 3+36*xi^2);

% Regularization for overlapping particles
regpoly = -1./(4*pi*r.^3) + 1/(4*pi)*(1-9*r/16+r.^3/32);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_1 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi) + (r < 2).*regpoly;


%% Field/Dipole coupling: rr component

% Polynomials multiplying the exponentials
exppolyp = 1./(512*pi^(3/2)*xi^5*r.^3).*(8*xi^4*r.^5 - 16*xi^4*r.^4 + 2*xi^2*(7-20*xi^2)*r.^3 - ...
    4*xi^2*(3-4*xi^2)*r.^2 - (3-12*xi^2+32*xi^4)*r - 2*(3+4*xi^2-32*xi^4));
exppolym = 1./(512*pi^(3/2)*xi^5*r.^3).*(8*xi^4*r.^5 + 16*xi^4*r.^4 + 2*xi^2*(7-20*xi^2)*r.^3 + ...
    4*xi^2*(3-4*xi^2)*r.^2 - (3-12*xi^2+32*xi^4)*r + 2*(3+4*xi^2-32*xi^4));
exppoly0 = 1./(256*pi^(3/2)*xi^5*r.^2).*(-8*xi^4*r.^4 - 2*xi^2*(7-36*xi^2)*r.^2 + 3-36*xi^2);

% Polynomials multiplying the error functions
erfpolyp = 1./(1024*pi*xi^6*r.^3).*(-16*xi^6*r.^6 - 36*xi^4*(1-4*xi^2)*r.^4 + 128*xi^6*r.^3 - 3+36*xi^2-256*xi^6);
erfpolym = 1./(1024*pi*xi^6*r.^3).*(-16*xi^6*r.^6 - 36*xi^4*(1-4*xi^2)*r.^4 - 128*xi^6*r.^3 - 3+36*xi^2-256*xi^6);
erfpoly0 = 1./(512*pi*xi^6*r.^3).*(16*xi^6*r.^6 + 36*xi^4*(1-4*xi^2)*r.^4 + 3-36*xi^2);

% Regularization for overlapping particles
regpoly = 1./(2*pi*r.^3) + 1/(4*pi)*(1-9*r./8+r.^3/8);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_2 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi) + (r < 2).*regpoly;

%% Field/Dipole Force: coefficient multiplying -(mi*mj)r and -( (mj*r)mi + (mi*r)mj - 2(mi*r)(mj*r)r )

% Polynomials multiplying the exponentials
exppolyp = 3./(1024*pi^(3/2)*xi^5*r.^4).*(4*xi^4*r.^5 - 8*xi^4*r.^4 + 4*xi^2*(1-2*xi^2)*r.^3 + 16*xi^4*r.^2 - (3-12*xi^2+32*xi^4)*r - 2*(3+4*xi^2-32*xi^4));
exppolym = 3./(1024*pi^(3/2)*xi^5*r.^4).*(4*xi^4*r.^5 + 8*xi^4*r.^4 + 4*xi^2*(1-2*xi^2)*r.^3 - 16*xi^4*r.^2 - (3-12*xi^2+32*xi^4)*r + 2*(3+4*xi^2-32*xi^4));
exppoly0 = 3./(512*pi^(3/2)*xi^5*r.^3).*(-4*xi^4*r.^4 - 4*xi^2*(1-6*xi^2)*r.^2 + 3-36*xi^2);

% Polynomials multiplying the error functions
erfpolyp = 3./(2048*pi*xi^6*r.^4).*(-8*xi^6*r.^6 - 12*xi^4*(1-4*xi^2)*r.^4 + 6*xi^2*(1-8*xi^2)*r.^2 - 3+36*xi^2-256*xi^6);
erfpolym = 3./(2048*pi*xi^6*r.^4).*(-8*xi^6*r.^6 - 12*xi^4*(1-4*xi^2)*r.^4 + 6*xi^2*(1-8*xi^2)*r.^2 - 3+36*xi^2-256*xi^6);
erfpoly0 = 3./(1024*pi*xi^6*r.^4).*(8*xi^6*r.^6 + 12*xi^4*(1-4*xi^2)*r.^4 - 6*xi^2*(1-8*xi^2)*r.^2 + 3-36*xi^2);

% Regularization for overlapping particles
regpoly = 3./(4*pi*r.^4) - 3/(64*pi)*(3-r.^2/2);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_force_1 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi) + (r < 2).*regpoly;

%% Field/Dipole Force from:  coefficient multiplying -(mi*r)(mj*r)r

% Polynomials multiplying the exponentials
exppolyp = 9./(1024*pi^(3/2)*xi^5*r.^4).*(4*xi^4*r.^5 - 8*xi^4*r.^4 + 8*xi^4*r.^3 + 8*xi^2*(1-2*xi^2)*r.^2 + (3-12*xi^2+32*xi^4)*r + 2*(3+4*xi^2-32*xi^4));
exppolym = 9./(1024*pi^(3/2)*xi^5*r.^4).*(4*xi^4*r.^5 + 8*xi^4*r.^4 + 8*xi^4*r.^3 - 8*xi^2*(1-2*xi^2)*r.^2 + (3-12*xi^2+32*xi^4)*r - 2*(3+4*xi^2-32*xi^4));
exppoly0 = 9./(512*pi^(3/2)*xi^5*r.^3).*(-4*xi^4*r.^4 + 8*xi^4*r.^2 - 3+36*xi^2);

% Polynomials multiplying the error functions
erfpolyp = 9./(2048*pi*xi^6*r.^4).*(-8*xi^6*r.^6 - 4*xi^4*(1-4*xi^2)*r.^4 - 2*xi^2*(1-8*xi^2)*r.^2 + 3-36*xi^2+256*xi^6);
erfpolym = 9./(2048*pi*xi^6*r.^4).*(-8*xi^6*r.^6 - 4*xi^4*(1-4*xi^2)*r.^4 - 2*xi^2*(1-8*xi^2)*r.^2 + 3-36*xi^2+256*xi^6);
erfpoly0 = 9./(1024*pi*xi^6*r.^4).*(8*xi^6*r.^6 + 4*xi^4*(1-4*xi^2)*r.^4 + 2*xi^2*(1-8*xi^2)*r.^2 - 3+36*xi^2);

% Regularization for overlapping particles
regpoly = -9./(4*pi*r.^4) - 9/(64*pi)*(1-r.^2/2);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_force_2 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi) + (r < 2).*regpoly;

%% Self terms

% Potential/charge
self = (1-exp(-4*xi^2))/(8*pi^(3/2)*xi) + erfc(2*xi)/(4*pi);
pot_charge = [self;pot_charge];

% Potential/dipole or field/charge
pot_dip = [0;pot_dip];

% Field/dipole
self = (-1+6*xi^2+(1-2*xi^2)*exp(-4*xi^2))/(16*pi^(3/2)*xi^3) + erfc(2*xi)/(4*pi);
field_dip_1 = [self;field_dip_1];
field_dip_2 = [self;field_dip_2];

% Field/dipole force
field_dip_force_1 = [0;field_dip_force_1];
field_dip_force_2 = [0;field_dip_force_2];

end