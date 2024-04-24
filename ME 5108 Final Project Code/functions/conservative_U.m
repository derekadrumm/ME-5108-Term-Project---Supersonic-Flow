%%%% Compute the conservative U vector from the primitaves %%%%%%%%%%%%%%%
% primitive variables: rho, U, V, Et

function [U1,U2,U3,U4] = conservative_U(rho,U,V,T,params)
paramcell = num2cell(params);
[cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl] = paramcell{:};

U1 = rho;

U2 = rho.*U;

U3 = rho.*V;

U4 = rho.*(cv.*T + (1/2)*(U.^2 + V.^2));

end