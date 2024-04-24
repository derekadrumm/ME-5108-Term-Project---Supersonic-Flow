%%%% Compute Cell Reynolds numbers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rex, Rey] = cell_reynolds(delx,dely,U,V,rho,mu)

Rex = (rho.*U*delx)./(mu);

Rey = (rho.*V*dely)./(mu);

end
