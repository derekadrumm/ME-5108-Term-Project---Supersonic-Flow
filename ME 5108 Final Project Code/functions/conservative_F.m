%%%% Compute the conservative E vector from the primitives %%%%%%%%%%%%%%%
% primitive variables: rho, U, V, P, Tauxy, Tauyy, Einternal, Qy

function [F1,F2,F3,F4] = conservative_F(rho,U,V,P,T,Tauxy,Tauyy,Qy,params)
paramcell = num2cell(params);
[cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl] = paramcell{:};

F1 = rho.*V;

F2 = rho.*U.*V - Tauxy;

F3 = rho.*(V.^2) + P - Tauyy;

Etotal = rho.*(cv.*T + (1/2)*(U.^2 + V.^2));
F4 = (Etotal + P).*V - U.*Tauxy - V.*Tauyy + Qy;

end