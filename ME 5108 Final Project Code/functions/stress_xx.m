%%%% Compute the Stress Tensor Tau %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute using a different finite difference (choice will depend on which
% step in MacCormacks method we are in!)
% variables needed: U, V, mu
% case variable - needed so that can change which finite difference to use
% case_var should be a STRING!

% to compute along boundaries, use forward difference for left and lower,
% and backward difference for right and upper!!

% case 1: 'Predict_E' (predictor step computation)
% case 2: 'Correct_E' (corrector step computation)

function Tauxx = stress_xx(delx,dely,U,V,mu,case_var)
Nx = size(U,2); Ny = size(U,1);
dudx = zeros(size(U)); dvdy = dudx;

switch case_var
    %%%% Backward in x, Central in y    
    case 'Predict_E'
        %%%% INTERIOR AND RIGHT BOUNDARY
        for i = 2:Ny-1
            for j = 2:Nx
           
                dudx(i,j) = (U(i,j)-U(i,j-1))/delx;
                dvdy(i,j) = (V(i+1,j)-V(i-1,j))/(2*dely);
                
            end
        end
        %%%% LEFT,BOTTOM,TOP BOUNDARIES
        dudx(:,1) = (U(:,2)-U(:,1))/delx;
        dvdy(1,:) = (V(2,:)-V(1,:))/dely;
        dvdy(end,:) = (V(end,:)-V(end-1,:))/dely;
        
        Tauxx = -(2/3)*mu.*(dudx + dvdy) + 2*mu.*dudx;
        
    %%%% Forward in x, Central in y        
    case 'Correct_E'
        %%%% INTERIOR AND LEFT BOUNDARY
        for i = 2:Ny-1
            for j = 1:Nx-1
           
                dudx(i,j) = (U(i,j+1)-U(i,j))/delx;
                dvdy(i,j) = (V(i+1,j)-V(i-1,j))/(2*dely);
                
            end
        end
        %%%% RIGHT,BOTTOM,TOP BOUNDARIES
        dudx(:,end) = (U(:,end)-U(:,end-1))/delx;
        dvdy(1,:) = (V(2,:)-V(1,:))/dely;
        dvdy(end,:) = (V(end,:)-V(end-1,:))/dely;
        
        Tauxx = -(2/3)*mu.*(dudx + dvdy) + 2*mu.*dudx;
        

%%%%%%%% Error case
    otherwise
        error('No finite difference scheme specified / incorrect "case_var"')
end
end
