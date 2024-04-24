%%%% Compute the Heat Flux Qx,Qy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute using a different finite difference (choice will depend on which
% step in MacCormacks method we are in!)
% variables needed: delx, dely, ktherm, T
% case variable - needed so that can change which finite difference to use
% case_var should be a STRING!

% to compute along boundaries, use forward difference for left and lower,
% and backward difference for right and upper!!

% case 1: FWDx_CENy
% case 2: CENx_FWDy
% case 3: BWDx_CENy
% case 4: CENx_BWDy

function [Qx,Qy] = heat_flux(delx,dely,ktherm,T,case_var)
Nx = size(T,2); Ny = size(T,1);
Qx = zeros(size(T)); Qy = Qx;

%%%% BOUNDARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qx(:,1) = -ktherm(:,1).*((T(:,2)-T(:,1))/delx);
Qx(:,end) = -ktherm(:,end).*((T(:,end)-T(:,end-1))/delx);

Qy(1,:) = -ktherm(1,:).*((T(2,:)-T(1,:))/dely);
Qy(end,:) = -ktherm(end,:).*((T(end,:)-T(end-1,:))/dely);

%%%% INTERIOR (CHANGES WITH PREDICTOR vs CORRECTOR STEP %%%%%%%%%%%%%%%%%%
switch case_var
%%%%%%%% Forward in x, Central in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'FWDx_CENy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                Qx(i,j) = -ktherm(i,j)*((T(i+1,j)-T(i,j))/(delx));
                Qy(i,j) = -ktherm(i,j)*((T(i,j+1)-T(i,j-1))/(2*dely));
                
            end
        end
             
%%%%%%%% Central in x, Forward in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    case 'CENx_FWDy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                Qx(i,j) = -ktherm(i,j)*((T(i+1,j)-T(i-1,j))/(2*delx));
                Qy(i,j) = -ktherm(i,j)*((T(i,j+1)-T(i,j))/(dely));
                
            end
        end

%%%%%%%% Backward in x, Central in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'BWDx_CENy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                Qx(i,j) = -ktherm(i,j)*((T(i,j)-T(i-1,j))/(delx));
                Qy(i,j) = -ktherm(i,j)*((T(i,j+1)-T(i,j-1))/(2*dely));
                
            end
        end
        
%%%%%%%% Central in x, Backward in y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'CENx_BWDy'
        for i = 2:Nx-1
            for j = 2:Ny-1
                
                Qx(i,j) = -ktherm(i,j)*((T(i+1,j)-T(i-1,j))/(2*delx));
                Qy(i,j) = -ktherm(i,j)*((T(i,j)-T(i,j-1))/(dely));
                
            end
        end

%%%%%%%% Error case
    otherwise
        error('No finite difference scheme specified / incorrect "case_var"')
end
end