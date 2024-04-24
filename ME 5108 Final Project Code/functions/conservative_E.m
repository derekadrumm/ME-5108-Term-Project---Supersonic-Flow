%%%% Compute the conservative E vector from the primitives %%%%%%%%%%%%%%%
% primitive variables: rho, U, V, P, Tauxx, Tauxy, Etotal, Qx

function [E1,E2,E3,E4] = conservative_E(rho,U,V,P,T,Tauxx,Tauxy,Qx,params)
paramcell = num2cell(params);
[cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl] = paramcell{:};

E1 = rho.*U;

E2 = rho.*(U.^2) + P - Tauxx;

E3 = rho.*U.*V - Tauxy;

Etotal = rho.*(cv.*T + (1/2)*(U.^2 + V.^2));
E4 = (Etotal+P).*U - U.*Tauxx - V.*Tauxy + Qx;

end