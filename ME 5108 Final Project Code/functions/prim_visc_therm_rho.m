%%%% Compute primitive viscosity thermal conductivity and density
% from the full T and P including boundary conditions
% primitive variables: mu, ktherm, rho
% parameter vector: params = [cv, cp, R, mu0, T0, prandtl]

function [mu,ktherm,rho] = prim_visc_therm_rho(P,T,params)
paramcell = num2cell(params);
[cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl] = paramcell{:};

rho = P./(R*T);

mu = (mu0*(T/Tfar).^(3/2)).*((Tfar+110)./(T+110));

ktherm = mu*cp/prandtl;
end