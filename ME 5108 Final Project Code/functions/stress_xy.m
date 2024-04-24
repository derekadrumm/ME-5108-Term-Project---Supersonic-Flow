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
% case 3: 'Predict_F' (predictor step computation)
% case 4: 'Correct_F' (corrector step computation)

function Tauxy = stress_xy(delx,dely,U,V,mu,case_var)
Nx = size(U,2); Ny = size(U,1);
dudy = zeros(size(U)); dvdx = dudy;

switch case_var
    %%%% Backward in x, Central in y    
    case 'Predict_E'
        %%%% INTERIOR AND RIGHT BOUNDARY
        for i = 2:Ny-1
            for j = 2:Nx
           
                dudy(i,j) = (U(i+1,j)-U(i-1,j))/(2*dely);
                dvdx(i,j) = (V(i,j)-V(i,j-1))/(delx);
                
            end
        end
        %%%% LEFT,BOTTOM,TOP BOUNDARIES
        dvdx(:,1) = (V(:,2)-V(:,1))/(delx);
        dudy(1,:) = (U(2,:)-U(1,:))/(dely);
        dudy(end,:) = (U(end,:)-U(end-1,:))/(dely);
        
        Tauxy = mu.*(dudy + dvdx);
        
        
    %%%% Forward in x, Central in y        
    case 'Correct_E'
        %%%% INTERIOR AND LEFT BOUNDARY
        for i = 2:Ny-1
            for j = 1:Nx-1
           
                dudy(i,j) = (U(i+1,j)-U(i-1,j))/(2*dely);
                dvdx(i,j) = (V(i,j+1)-V(i,j))/(delx);
                
            end
        end
        %%%% RIGHT,BOTTOM,TOP BOUNDARIES
        dvdx(:,end) = (V(:,end)-V(:,end-1))/(delx);
        dudy(1,:) = (U(2,:)-U(1,:))/(dely);
        dudy(end,:) = (U(end,:)-U(end-1,:))/(dely);
        
        Tauxy = mu.*(dudy + dvdx);
        
        
    %%%% Central in x, Backward in y    
    case 'Predict_F'
        %%%% INTERIOR AND BOTTOM BOUNDARY
        for i = 2:Ny
            for j = 2:Nx-1
           
                dudy(i,j) = (U(i,j)-U(i-1,j))/(dely);
                dvdx(i,j) = (V(i,j+1)-V(i,j-1))/(2*delx);
                
            end
        end
        %%%% LEFT,RIGHT,TOP BOUNDARIES
        dvdx(:,1) = (V(:,2)-V(:,1))/(delx);
        dvdx(:,end) = (V(:,end)-V(:,end-1))/(delx);
        dudy(1,:) = (U(2,:)-U(1,:))/(dely);
        
        Tauxy = mu.*(dudy + dvdx);
        
        
    %%%% Central in x, Forward in y        
    case 'Correct_F'
        %%%% INTERIOR AND TOP BOUNDARY
        for i = 1:Ny-1
            for j = 2:Nx-1
           
                dudy(i,j) = (U(i+1,j)-U(i,j))/(dely);
                dvdx(i,j) = (V(i,j+1)-V(i,j-1))/(2*delx);
                
            end
        end
        %%%% LEFT,RIGHT,BOTTOM BOUNDARIES
        dvdx(:,1) = (V(:,2)-V(:,1))/(delx);
        dvdx(:,end) = (V(:,end)-V(:,end-1))/(delx);
        dudy(end,:) = (U(end,:)-U(end-1,:))/(dely);
        
        Tauxy = mu.*(dudy + dvdx);

%%%%%%%% Error case
    otherwise
        error('No finite difference scheme specified / incorrect "case_var"')
end
end