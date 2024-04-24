%%%% Compute the Stress Tensor Tau %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute using a different finite difference (choice will depend on which
% step in MacCormacks method we are in!)
% variables needed: U, V, mu
% case variable - needed so that can change which finite difference to use
% case_var should be a STRING!

% to compute along boundaries, use forward difference for left and lower,
% and backward difference for right and upper!!

% case 1: FWDx_CENy
% case 2: CENx_FWDy
% case 3: BWDx_CENy
% case 4: CENx_BWDy

function [Tauxx, Tauxy, Tauyy] = stress_tensor(delx,dely,U,V,mu,case_var)
Nx = size(U,2); Ny = size(U,1);
Tauxx = zeros(size(U)); Tauxy = Tauxx; Tauyy = Tauxx;

%%%% BOUNDARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tauxx(:,1) = -2/3 * mu(:,1)*((U(:,2)-U(:,1))/delx + (V(:,2)-V(:,1))/(2*dely)) + ...
    ;

%%%% INTERIOR (CHANGES WITH PREDICTOR vs CORRECTOR STEP %%%%%%%%%%%%%%%%%%
switch case_var
%%%%%%%% Forward in x, Central in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'FWDx_CENy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                divUV = -2/3 * mu(i,j)*((U(i+1,j)-U(i,j))/delx + (V(i,j+1)-V(i,j-1))/(2*dely));
                
                Tauxx(i,j) = divUV + 2*mu(i,j)*((U(i+1,j)-U(i,j))/delx);
                
                Tauxy(i,j) = mu(i,j)*((U(i,j+1)-U(i,j-1))/(2*dely) + (V(i+1,j)-V(i,j))/delx);
                
                Tauyy(i,j) = divUV + 2*mu(i,j)*((V(i,j+1)-V(i,j-1))/(2*dely));
                
            end
        end
             
%%%%%%%% Central in x, Forward in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    case 'CENx_FWDy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                divUV = -2/3 * mu(i,j)*((U(i+1,j)-U(i-1,j))/(2*delx) + (V(i,j+1)-V(i,j))/(dely));
                
                Tauxx(i,j) = divUV + 2*mu(i,j)*((U(i+1,j)-U(i-1,j))/(2*delx));
                
                Tauxy(i,j) = mu(i,j)*((U(i,j+1)-U(i,j))/(dely) + (V(i+1,j)-V(i-1,j))/(2*delx));
                
                Tauyy(i,j) = divUV + 2*mu(i,j)*((V(i,j+1)-V(i,j))/(dely));
                
            end
        end

%%%%%%%% Backward in x, Central in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'BWDx_CENy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                divUV = -2/3 * mu(i,j)*((U(i,j)-U(i-1,j))/(delx) + (V(i,j+1)-V(i,j-1))/(2*dely));
                
                Tauxx(i,j) = divUV + 2*mu(i,j)*((U(i,j)-U(i-1,j))/(delx));
                
                Tauxy(i,j) = mu(i,j)*((U(i,j+1)-U(i,j-1))/(2*dely) + (V(i,j)-V(i-1,j))/(delx));
                
                Tauyy(i,j) = divUV + 2*mu(i,j)*((V(i,j+1)-V(i,j-1))/(2*dely));
                
            end
        end
        
%%%%%%%% Central in x, Backward in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'CENx_BWDy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                divUV = -2/3 * mu(i,j)*((U(i+1,j)-U(i-1,j))/(2*delx) + (V(i,j)-V(i,j-1))/(dely));
                
                Tauxx(i,j) = divUV + 2*mu(i,j)*((U(i+1,j)-U(i-1,j))/(2*delx));
                
                Tauxy(i,j) = mu(i,j)*((U(i,j)-U(i,j-1))/(dely) + (V(i+1,j)-V(i-1,j))/(2*delx));
                
                Tauyy(i,j) = divUV + 2*mu(i,j)*((V(i,j)-V(i,j-1))/(dely));
                
            end
        end

%%%%%%%% Error case
    otherwise
        error('No finite difference scheme specified / incorrect "case_var"')
end
end
