%%%% Compute the Stress Tensor Tau %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute using a different finite difference (choice will depend on which
% step in MacCormacks method we are in!)
% variables needed: U, V, mu
% case variable - needed so that can change which finite difference to use
% case_var should be a STRING!

% to compute along boundaries, use forward difference for left and lower,
% and backward difference for right and upper!!

% case 1: 'Predict_F' (predictor step computation)
% case 2: 'Correct_F' (corrector step computation)

function Tauyy = stress_yy(delx,dely,U,V,mu,case_var)
Nx = size(U,2); Ny = size(U,1);
dudx = zeros(size(U)); dvdy = dudx;

switch case_var
%%%% Central in x, Backward in y    
    case 'Predict_F'
        %%%% INTERIOR AND BOTTOM BOUNDARY
        for i = 2:Ny
            for j = 2:Nx-1
           
                dudx(i,j) = (U(i,j+1)-U(i,j-1))/(2*delx);
                dvdy(i,j) = (V(i,j)-V(i-1,j))/(dely);
                
            end
        end
        %%%% LEFT,RIGHT,TOP BOUNDARIES
        dudx(:,1) = (U(:,2)-U(:,1))/delx;
        dudx(:,end) = (U(:,end)-U(:,end-1))/delx;
        dvdy(1,:) = (V(2,:)-V(1,:))/dely;
        
        Tauyy = -(2/3)*mu.*(dudx + dvdy) + 2*mu.*dvdy;
        
%%%% Central in x, Forward in y        
    case 'Correct_F'
        %%%% INTERIOR AND TOP BOUNDARY
        for i = 1:Ny-1
            for j = 2:Nx-1
           
                dudx(i,j) = (U(i,j+1)-U(i,j-1))/(2*delx);
                dvdy(i,j) = (V(i+1,j)-V(i,j))/(dely);
                
            end
        end
        %%%% LEFT,RIGHT,BOTTOM BOUNDARIES
        dudx(:,1) = (U(:,2)-U(:,1))/delx;
        dudx(:,end) = (U(:,end)-U(:,end-1))/delx;
        dvdy(end,:) = (V(end,:)-V(end-1,:))/dely;
        
        Tauyy = -(2/3)*mu.*(dudx + dvdy) + 2*mu.*dvdy;
        

%%%%%%%% Error case
    otherwise
        error('No finite difference scheme specified / incorrect "case_var"')
end
end