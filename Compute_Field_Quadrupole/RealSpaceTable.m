function [pot_charge,pot_dip,pot_quad,field_dip_1,field_dip_2,field_dip_force_1,...
         field_dip_force_2,field_quad_1,field_quad_2,field_quad_3,grad_quad_1, ...
         grad_quad_2, grad_quad_3,grad_quad_4] = RealSpaceTable(r,xi)

% Tabulate the real space contributions as a function of particle separation 
%
% INPUTS
% r = (N-by-1) separation values to tabulate
% xi = (scalar) Ewald splitting parameter
%
% OUTPUTS
% pot_charge = (N-by-1) potential/charge coupling
% pot_dip = (N-by-1) potential/dipole or field/charge coupling
% pot_quad = (N-by-1) potential/quadrupole or field gradient/charge coupling
% field_dip_1 = (N-by-1) field/dipole coupling: I-rr component
% field_dip_2 = (N-by-1) field/dipole coupling: rr component
% field_dip_force_1 = (N-by-1) field/dipole force: -(Si*Sj)r and -((Sj*r)Si + (Si*r)Sj - 2(Si*r)(Sj*r)r) components
% field_dip_force_2 = (N-by-1) field/dipole force: (Si*r)(Sj*r)r component
% field_quad_1 = (N-by-1) field/quadrupole or field gradient/dipole coupling: rrr*S component
% field_quad_2 = (N-by-1) field/quadrupole or field gradient/dipole coupling: Ir*S+Sr+rS component
% field_quad_3 = (N-by-1) field/quadrupole or field gradient/dipole coupling: Ir*S component
% grad_quad_1 = (N-by-1) field gradient/quadrupole coupling: Irr:Q component
% grad_quad_2 = (N-by-1) field gradient/quadrupole coupling: Q component
% grad_quad_3 = (N-by-1) field gradient/quadrupole coupling: Irr:Q+2Q*rr+2rr*Q component
% grad_quad_4 = (N-by-1) field gradient/quadrupole coupling: rrrr:Q component

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
%regpoly = -1./(4*pi*r) + (4-r)./(16*pi);

% Combine the polynomial coefficients, exponentials, and error functions
pot_charge = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

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
%regpoly = -1./(4*pi*r.^2) + r/(8*pi).*(1-3/8*r);

% Combine the polynomial coefficients, exponentials, and error functions
pot_dip = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Potential/Quadrupole or Field Gradient/Charge Coupling

% Polynomials multiplying the exponetials
exppolyp = 3./(1024*pi^(3/2)*xi^5*r.^3).*(20*xi^4*r.^5 + 24*xi^4*r.^4 + 2*xi^2*(10-4*xi^2)*r.^3 ...
    - 16*xi^2*(2-xi^2)*r.^2 - (15-12*xi^2+32*xi^4)*r + 2*(9-4*xi^2+32*xi^4));
exppolym = 3./(1024*pi^(3/2)*xi^5*r.^3).*(20*xi^4*r.^5 - 24*xi^4*r.^4 + 2*xi^2*(10-4*xi^2)*r.^3 ...
    + 16*xi^2*(2-xi^2)*r.^2 - (15-12*xi^2+32*xi^4)*r - 2*(9-4*xi^2+32*xi^4));
exppoly0 = 15./(512*pi^(3/2)*xi^5*r.^2).*(4*xi^4*r.^4 + 4*xi^2*(1+2*xi^2)*r.^2 - 3*(1+4*xi^2));

% Polynomials multiplying the error functions
erfpolyp = 3./(2048*pi*xi^6*r.^3).*(40*xi^6*r.^6 + 128*xi^6*r.^5 + 20*xi^4*(3+4*xi^2)*r.^4 ...
    - 10*xi^2*(3+8*xi^2)*r.^2 + 15+60*xi^2+256*xi^6);
erfpolym = 3./(2048*pi*xi^6*r.^3).*(40*xi^6*r.^6 - 128*xi^6*r.^5 + 20*xi^4*(3+4*xi^2)*r.^4 ...
    - 10*xi^2*(3+8*xi^2)*r.^2 + 15+60*xi^2+256*xi^6);
erfpoly0 = 15./(1024*pi*xi^6*r.^3).*(-8*xi^6*r.^6 - 4*xi^4*(3+4*xi^2)*r.^4 - 2*xi^2*(3+8*xi^2)*r.^2 - 3*(1+4*xi^2));

% Regularization for overlapping particles
%regpoly = 15./(1024*pi*xi^6*r.^3).*(12*xi^4*r.^4 - 2*xi^2*(3+8*xi^2)*r.^2 + 3+12*xi^2); % this is wrong; shouldn't have any xi in it

% Combine the polynomial coefficients, exponentials, and error functions
pot_quad = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

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
%regpoly = -1./(4*pi*r.^3) + 1/(4*pi)*(1-9*r/16+r.^3/32);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_1 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);


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
%regpoly = 1./(2*pi*r.^3) + 1/(4*pi)*(1-9*r./8+r.^3/8);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_2 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

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
%regpoly = 3./(4*pi*r.^4) - 3/(64*pi)*(3-r.^2/2);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_force_1 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

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
%regpoly = -9./(4*pi*r.^4) - 9/(64*pi)*(1-r.^2/2);

% Combine the polynomial coefficients, exponentials, and error functions
field_dip_force_2 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Field/Quadrupole or Field Gradient/Dipole Coupling: rrr*S component

% Polynomials multiplying the exponentials
exppolyp = 15./(16384*pi^(3/2)*xi^7*r.^4).*( -24*xi^6*r.^7 + 48*xi^6*r.^6 - 4*xi^4*(9-8*xi^2).*r.^5 + 8*xi^4*(3-8*xi^2).*r.^4 ...
    + 2*xi^2*(21-80*xi^2+64*xi^4)*r.^3 - 4*xi^2*(3-8*xi^2)^2*r.^2 - (45-120*xi^2+192*xi^4-512*xi^6)*r - 2*(45+24*xi^2-64*xi^4+512*xi^6) );
exppolym = 15./(16384*pi^(3/2)*xi^7*r.^4).*( -24*xi^6*r.^7 - 48*xi^6*r.^6 - 4*xi^4*(9-8*xi^2).*r.^5 - 8*xi^4*(3-8*xi^2).*r.^4 ...
    + 2*xi^2*(21-80*xi^2+64*xi^4)*r.^3 + 4*xi^2*(3-8*xi^2)^2*r.^2 - (45-120*xi^2+192*xi^4-512*xi^6)*r + 2*(45+24*xi^2-64*xi^4+512*xi^6) );
exppoly0 = 15./(8192*pi^(3/2)*xi^7*r.^3).*( 24*xi^6*r.^6 + 4*xi^4*(9-32*xi^2)*r.^4 - 2*xi^2*(21-128*xi^2)*r.^2 + 45-480*xi^2 );

% Polynomials multiplying the error functions
erfpolyp = 15./(32768*pi*xi^8*r.^4).*( 48*xi^8*r.^8 + 32*xi^6*(3-8*xi^2)*r.^6 - 24*xi^4*(3-16*xi^2)*r.^4 ...
    + 72*xi^2*(1-8*xi^2)*r.^2 -45+480*xi^2+4096*xi^8 );
erfpolym = 15./(32768*pi*xi^8*r.^4).*( 48*xi^8*r.^8 + 32*xi^6*(3-8*xi^2)*r.^6 - 24*xi^4*(3-16*xi^2)*r.^4 ...
    + 72*xi^2*(1-8*xi^2)*r.^2 -45+480*xi^2+4096*xi^8 );
erfpoly0 = 15./(16384*pi*xi^8*r.^4).*( -48*xi^8*r.^8 - 32*xi^6*(3-8*xi^2)*r.^6 + 24*xi^4*(3-16*xi^2)*r.^4 ...
    - 72*xi^2*(1-8*xi^2)*r.^2 +45-480*xi^2 );

% Regularization for overlapping particles
%regpoly = -15./(4*pi*r.^4) + 15*r.^2/(64*pi).*(1-3*r.^2/16);

% Combine the polynomial coefficients, exponentials, and error functions
field_quad_1 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Field/Quadrupole or Field Gradient/Dipole Coupling: Ir*S+Sr+rS component

% Polynomials multiplying the exponentials
exppolyp = 3./(16384*pi^(3/2)*xi^7*r.^4).*( -40*xi^6*r.^7 + 80*xi^6*r.^6 - 20*xi^4*(11-24*xi^2)*r.^5 + 8*xi^4*(45+8*xi^2)*r.^4 ...
    - 2*xi^2*(45-80*xi^2+64*xi^4)*r.^3 - 4*xi^2*(15+48*xi^2-64*xi^4)*r.^2 + (45-120*xi^2+192*xi^4-512*xi^6)*r + 2*(45+24*xi^2-64*xi^4+512*xi^6) );
exppolym = 3./(16384*pi^(3/2)*xi^7*r.^4).*( -40*xi^6*r.^7 - 80*xi^6*r.^6 - 20*xi^4*(11-24*xi^2)*r.^5 - 8*xi^4*(45+8*xi^2)*r.^4 ...
    - 2*xi^2*(45-80*xi^2+64*xi^4)*r.^3 + 4*xi^2*(15+48*xi^2-64*xi^4)*r.^2 + (45-120*xi^2+192*xi^4-512*xi^6)*r - 2*(45+24*xi^2-64*xi^4+512*xi^6) );
exppoly0 = 15./(8192*pi^(3/2)*xi^7*r.^3).*( 8*xi^6*r.^6 + 4*xi^4*(11-32*xi^2)*r.^4 + 2*xi^2*(9-64*xi^2)*r.^2 -9+96*xi^2 );

% Polynomials multiplying the error functions
erfpolyp = 3./(32768*pi*xi^8*r.^4).*( 80*xi^8*r.^8 + 160*xi^6*(3-8*xi^2)*r.^6 - 2048*xi^8*r.^5 + 120*xi^4*(3-16*xi^2).*r.^4 ...
    - 120*xi^2*(1-8*xi^2)*r.^2 + 45-480*xi^2-4096*xi^8 );
erfpolym = 3./(32768*pi*xi^8*r.^4).*( 80*xi^8*r.^8 + 160*xi^6*(3-8*xi^2)*r.^6 + 2048*xi^8*r.^5 + 120*xi^4*(3-16*xi^2).*r.^4 ...
    - 120*xi^2*(1-8*xi^2)*r.^2 + 45-480*xi^2-4096*xi^8 );
erfpoly0 = 15./(16384*pi*xi^8*r.^4).*( -16*xi^8*r.^8 - 32*xi^6*(3-8*xi^2)*r.^6 - 24*xi^4*(3-16*xi^2)*r.^4 + 24*xi^2*(1-8*xi^2)*r.^2 -9+96*xi^2);

% Regularization for overlapping particles
%regpoly = 3./(4*pi*r.^4) - 3*r/(8*pi).*(1-5*r/8+5*r.^3/128);

% Combine the polynomial coefficients, exponentials, and error functions
field_quad_2 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Field/Quadrupole or Field Gradient/Dipole Coupling: Ir*S component

% Polynomials multiplying the exponentials
exppolyp = 5./(1024*pi^(3/2)*xi^5*r.^2).*( 4*xi^4*r.^5 - 8*xi^4*r.^4 + 16*xi^2*(1-2*xi^2)*r.^3 - 24*xi^2*r.^2 + 3*r + 6);
exppolym = 5./(1024*pi^(3/2)*xi^5*r.^2).*( 4*xi^4*r.^5 + 8*xi^4*r.^4 + 16*xi^2*(1-2*xi^2)*r.^3 + 24*xi^2*r.^2 + 3*r - 6);
exppoly0 = 5./(512*pi^(3/2)*xi^5*r).*( -4*xi^4*r.^4 - 16*xi^2*(1-3*xi^2)*r.^2 -3+24*xi^2 );

% Polynomials multiplying the error functions
erfpolyp = 5./(2048*pi*xi^6*r.^2).*( -8*xi^6*r.^6 - 12*xi^4*(3-8*xi^2)*r.^4 + 128*xi^6*r.^3 - 6*xi^2*(3-16*xi^2)*r.^2 + 3-24*xi^2 );
erfpolym = 5./(2048*pi*xi^6*r.^2).*( -8*xi^6*r.^6 - 12*xi^4*(3-8*xi^2)*r.^4 - 128*xi^6*r.^3 - 6*xi^2*(3-16*xi^2)*r.^2 + 3-24*xi^2 );
erfpoly0 = 5./(1024*pi*xi^6*r.^2).*( 8*xi^6*r.^6 + 12*xi^4*(3-8*xi^2)*r.^4 + 6*xi^2*(3-16*xi^2)*r.^2 -3+24*xi^2 );

% Regularization for overlapping particles
%regpoly = 5*r/(8*pi).*(1-3*r/4+r.^3/16);

% Combine the polynomial coefficients, exponentials, and error functions
field_quad_3 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Field Gradient/Quadrupole Coupling: Irr:Q component

% Polynomials multiplying the exponentials
exppolyp = 75./(16384*pi^(3/2)*xi^7*r.^3).*( 8*xi^6*r.^7 - 16*xi^6*r.^6 + (44*xi^4-32*xi^6)*r.^5 - (72*xi^4-64*xi^6)*r.^4 ...
    + (18*xi^2+32*xi^4)*r.^3 + (12*xi^2-64*xi^4)*r.^2 - (9+24*xi^2)*r - (18-48*xi^2) );
exppolym = 75./(16384*pi^(3/2)*xi^7*r.^3).*( 8*xi^6*r.^7 + 16*xi^6*r.^6 + (44*xi^4-32*xi^6)*r.^5 + (72*xi^4-64*xi^6)*r.^4 ...
    + (18*xi^2+32*xi^4)*r.^3 - (12*xi^2-64*xi^4)*r.^2 - (9+24*xi^2)*r + (18-48*xi^2) );
exppoly0 = -75./(8192*pi^(3/2)*xi^7*r.^2).*( 8*xi^6*r.^6 + (44*xi^4-64*xi^6)*r.^4 + (18*xi^2-64*xi^4+128*xi^6)*r.^2 -9+48*xi^2-192*xi^4 );

% Polynomials multiplying the error functions
erfpolyp = -75./(32768*pi*xi^8*r.^3).*( 16*xi^8*r.^8 + (96*xi^6-128*xi^8)*r.^6 + (72*xi^4-192*xi^6+256*xi^8)*r.^4 ... 
    - (24*xi^2-96*xi^4+256*xi^6)*r.^2 +9-48*xi^2+192*xi^4 );
erfpolym = -75./(32768*pi*xi^8*r.^3).*( 16*xi^8*r.^8 + (96*xi^6-128*xi^8)*r.^6 + (72*xi^4-192*xi^6+256*xi^8)*r.^4 ... 
    - (24*xi^2-96*xi^4+256*xi^6)*r.^2 +9-48*xi^2+192*xi^4 );
erfpoly0 = 75./(16384*pi*xi^8*r.^3).*( 16*xi^8*r.^8 + (96*xi^6-128*xi^8)*r.^6 + (72*xi^4-192*xi^6+256*xi^8)*r.^4 ... 
    - (24*xi^2-96*xi^4+256*xi^6)*r.^2 +9-48*xi^2+192*xi^4 );

% Regularization for overlapping particles
%regpoly = 75*r.*(4-r.^2).^2/(1024*pi);

% Combine the polynomial coefficients, exponentials, and error functions
grad_quad_1 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Field Gradient/Quadrupole Coupling: Q component

% Polynomials multiplying the exponentials
exppolyp = 3./(32768*pi^(3/2)*xi^9*r.^5).*( -48*xi^8*r.^9 + 96*xi^8*r.^8 - (576*xi^6-608*xi^8)*r.^7 + (1056*xi^6-1216*xi^8)*r.^6 ...
    - (1536*xi^4-2576*xi^6+3968*xi^8)*r.^5 + (2160*xi^4-4320*xi^6-256*xi^8)*r.^4 - (360*xi^2+360*xi^4+640*xi^6-512*xi^8)*r.^3 ...
    - (360*xi^2-1680*xi^4-768*xi^6+1024*xi^8)*r.^2 + (135+180*xi^2+480*xi^4-768*xi^6+2048*xi^8)*r + (270-1080*xi^2-192*xi^4+512*xi^6-4096*xi^8) );
exppolym = 3./(32768*pi^(3/2)*xi^9*r.^5).*( -48*xi^8*r.^9 - 96*xi^8*r.^8 - (576*xi^6-608*xi^8)*r.^7 - (1056*xi^6-1216*xi^8)*r.^6 ...
    - (1536*xi^4-2576*xi^6+3968*xi^8)*r.^5 - (2160*xi^4-4320*xi^6-256*xi^8)*r.^4 - (360*xi^2+360*xi^4+640*xi^6-512*xi^8)*r.^3 ...
    + (360*xi^2-1680*xi^4-768*xi^6+1024*xi^8)*r.^2 + (135+180*xi^2+480*xi^4-768*xi^6+2048*xi^8)*r - (270-1080*xi^2-192*xi^4+512*xi^6-4096*xi^8) ); 
exppoly0 = 1./(16384*pi^(3/2)*xi^9*r.^4).*( 144*xi^8*r.^8 + (1728*xi^6-2400*xi^8)*r.^6 + (4608*xi^4-13200*xi^6+19200*xi^8)*r.^4 ...
    + (1080*xi^2-5400*xi^4+19200*xi^6)*r.^2 -405+2700*xi^2-14400*xi^4 );

% Polynomials multiplying the error functions
erfpolyp = 1./(65536*pi*xi^10*r.^5).*( 288*xi^10*r.^10 + (3600*xi^8-4800*xi^10)*r.^8 + (10800*xi^6-28800*xi^8+38400*xi^10)*r.^6 ...
    + 49152*xi^10*r.^5 + (5400*xi^4-21600*xi^6+57600*xi^8)*r.^4 - (1350*xi^2-7200*xi^4+28800*xi^6)*r.^2 +405-2700*xi^2+14400*xi^4+49152*xi^10 );
erfpolym = 1./(65536*pi*xi^10*r.^5).*( 288*xi^10*r.^10 + (3600*xi^8-4800*xi^10)*r.^8 + (10800*xi^6-28800*xi^8+38400*xi^10)*r.^6 ...
    - 49152*xi^10*r.^5 + (5400*xi^4-21600*xi^6+57600*xi^8)*r.^4 - (1350*xi^2-7200*xi^4+28800*xi^6)*r.^2 +405-2700*xi^2+14400*xi^4+49152*xi^10 );
erfpoly0 = -1./(32768*pi*xi^10*r.^5).*( 288*xi^10*r.^10 + (3600*xi^8-4800*xi^10)*r.^8 + (10800*xi^6-28800*xi^8+38400*xi^10)*r.^6 ...
    + (5400*xi^4-21600*xi^6+57600*xi^8)*r.^4 - (1350*xi^2-7200*xi^4+28800*xi^6)*r.^2 +405-2700*xi^2+14400*xi^4 );

% Regularization for overlapping particles
%regpoly = -3/(2*pi)*(1./r.^5-1+25*r/32-25*r.^3/256+3*r.^5/512);

% Combine the polynomial coefficients, exponentials, and error functions
grad_quad_2 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Field Gradient/Quadrupole Coupling: Irr:Q+2Q*rr+2rr*Q component

% Polynomials multiplying the exponentials
exppolyp = -15./(65536*pi^(3/2)*xi^9*r.^5).*( 48*xi^8*r.^9 - 96*xi^8*r.^8 + (336*xi^6-288*xi^8)*r.^7 - 576*(xi^6-xi^8)*r.^6 ...
    + (216*xi^4+144*xi^6+128*xi^8)*r.^5 - (480*xi^6+256*xi^8)*r.^4 - (180*xi^2-120*xi^4+640*xi^6-512*xi^8)*r.^3 ...
    + (720*xi^4+768*xi^6-1024*xi^8)*r.^2 + (135+180*xi^2+480*xi^4-768*xi^6+2048*xi^8)*r - (-270+1080*xi^2+192*xi^4-512*xi^6+4096*xi^8) );
exppolym = -15./(65536*pi^(3/2)*xi^9*r.^5).*( 48*xi^8*r.^9 + 96*xi^8*r.^8 + (336*xi^6-288*xi^8)*r.^7 + 576*(xi^6-xi^8)*r.^6 ...
    + (216*xi^4+144*xi^6+128*xi^8)*r.^5 + (480*xi^6+256*xi^8)*r.^4 - (180*xi^2-120*xi^4+640*xi^6-512*xi^8)*r.^3 ...
    - (720*xi^4+768*xi^6-1024*xi^8)*r.^2 + (135+180*xi^2+480*xi^4-768*xi^6+2048*xi^8)*r + (-270+1080*xi^2+192*xi^4-512*xi^6+4096*xi^8) );
exppoly0 = 15./(32768*pi^(3/2)*xi^9*r.^4).*( 48*xi^8*r.^8 + (336*xi^6-480*xi^8)*r.^6 + (216*xi^4-720*xi^6+1280*xi^8)*r.^4 ...
    - (180*xi^2-840*xi^4+2560*xi^6)*r.^2 +135-900*xi^2+4800*xi^4 );

% Polynomials multiplying the error functions
erfpolyp = -15./(131072*pi*xi^10*r.^5).*( -96*xi^10*r.^10 + (-720*xi^8+960*xi^10)*r.^8 - (720*xi^6-1920*xi^8+2560*xi^10)*r.^6 ...
    + (360*xi^4-1440*xi^6+3840*xi^8)*r.^4 - (270*xi^2-1440*xi^4+5760*xi^6)*r.^2 +135-900*xi^2+4800*xi^4+16384*xi^10 );
erfpolym = -15./(131072*pi*xi^10*r.^5).*( -96*xi^10*r.^10 + (-720*xi^8+960*xi^10)*r.^8 - (720*xi^6-1920*xi^8+2560*xi^10)*r.^6 ...
    + (360*xi^4-1440*xi^6+3840*xi^8)*r.^4 - (270*xi^2-1440*xi^4+5760*xi^6)*r.^2 +135-900*xi^2+4800*xi^4+16384*xi^10 );
erfpoly0 = -15./(65536*pi*xi^10*r.^5).*( 96*xi^10*r.^10 + (720*xi^8-960*xi^10)*r.^8 + (720*xi^6-1920*xi^8+2560*xi^10)*r.^6 ...
    - (360*xi^4-1440*xi^6+3840*xi^8)*r.^4 + (270*xi^2-1440*xi^4+5760*xi^6)*r.^2 -135+900*xi^2-4800*xi^4 );

% Regularization for overlapping particles
%regpoly = -15*(r.^2-4).^3.*(3*r.^4+6*r.^2+8)./(2048*pi*r.^5);

% Combine the polynomial coefficients, exponentials, and error functions
grad_quad_3 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

%% Field Gradient/Quadrupole Coupling: rrrr:Q component

% Polynomials multiplying the exponentials
exppolyp = -15./(65536*pi^(3/2)*xi^9*r.^5).*( 144*xi^8*r.^9 - 288*xi^8*r.^8 + (288*xi^6+96*xi^8)*r.^7 - (288*xi^6+192*xi^8)*r.^6 ...
    - (432*xi^4-912*xi^6+896*xi^8)*r.^5 + (720*xi^4-480*xi^6+1792*xi^8)*r.^4 + (720*xi^2-2280*xi^4+4480*xi^6-3584*xi^8)*r.^3 ...
    - (1080*xi^2+2160*xi^4+5376*xi^6-7168*xi^8)*r.^2 - (945+1260*xi^2+3360*xi^4-5376*xi^6+14336*xi^8)*r ...
    - (1890-7560*xi^2-1344*xi^4+3584*xi^6-28672*xi^8) );
exppolym = -15./(65536*pi^(3/2)*xi^9*r.^5).*( 144*xi^8*r.^9 + 288*xi^8*r.^8 + (288*xi^6+96*xi^8)*r.^7 + (288*xi^6+192*xi^8)*r.^6 ...
    - (432*xi^4-912*xi^6+896*xi^8)*r.^5 - (720*xi^4-480*xi^6+1792*xi^8)*r.^4 + (720*xi^2-2280*xi^4+4480*xi^6-3584*xi^8)*r.^3 ...
    + (1080*xi^2+2160*xi^4+5376*xi^6-7168*xi^8)*r.^2 - (945+1260*xi^2+3360*xi^4-5376*xi^6+14336*xi^8)*r ...
    + (1890-7560*xi^2-1344*xi^4+3584*xi^6-28672*xi^8) );
exppoly0 = 15./(32768*pi^(3/2)*xi^9*r.^4).*( 144*xi^8*r.^8 + (288*xi^6-480*xi^8)*r.^6 - (432*xi^4-1200*xi^6+1280*xi^8)*r.^4 ...
    + (720*xi^2-3000*xi^4+6400*xi^6)*r.^2 -945+6300*xi^2-33600*xi^4 );

% Polynomials multiplying the error functions
erfpolyp = 15./(131072*pi*xi^10*r.^5).*( 288*xi^10*r.^10 + (720*xi^8-960*xi^10)*r.^8 - (720*xi^6-1920*xi^8+2560*xi^10)*r.^6 ... 
    + (1080*xi^4-4320*xi^6+11520*xi^8)*r.^4 - (1350*xi^2-7200*xi^4+28800*xi^6)*r.^2 +945-6300*xi^2+33600*xi^4+114688*xi^10 ); 
erfpolym = 15./(131072*pi*xi^10*r.^5).*( 288*xi^10*r.^10 + (720*xi^8-960*xi^10)*r.^8 - (720*xi^6-1920*xi^8+2560*xi^10)*r.^6 ... 
    + (1080*xi^4-4320*xi^6+11520*xi^8)*r.^4 - (1350*xi^2-7200*xi^4+28800*xi^6)*r.^2 +945-6300*xi^2+33600*xi^4+114688*xi^10 );
erfpoly0 = -15./(65536*pi*xi^10*r.^5).*( 288*xi^10*r.^10 + (720*xi^8-960*xi^10)*r.^8 - (720*xi^6-1920*xi^8+2560*xi^10)*r.^6 ... 
    + (1080*xi^4-4320*xi^6+11520*xi^8)*r.^4 - (1350*xi^2-7200*xi^4+28800*xi^6)*r.^2 +945-6300*xi^2+33600*xi^4 );

% Regularization for overlapping particles
%regpoly = -15*(3584 - 80*r.^6 - 30*r.^8 + 9*r.^10)./(2048*pi*r.^5);

% Combine the polynomial coefficients, exponentials, and error functions
grad_quad_4 = exppolyp.*exp(-(r+2).^2*xi^2) + exppolym.*exp(-(r-2).^2*xi^2) + exppoly0.*exp(-r.^2*xi^2) + ...
    erfpolyp.*erfc((r+2)*xi) + erfpolym.*erfc((r-2)*xi) + erfpoly0.*erfc(r*xi);

end