%%%% Compute the primitive variables from the conservatives %%%%%%%%%%%%%%%
% primitive variables: U, V, P, T
% parameter vector: params = [cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl]

function [U,V,P,T] = primitive_UVPT(U1,U2,U3,U4,params)
paramcell = num2cell(params);
[cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl] = paramcell{:};



U = U2./U1;

V = U3./U1;

Einternal = U4./U1 - (1/2)*(U.^2 + V.^2);

T = Einternal/cv;

P = R*U1.*T;
end