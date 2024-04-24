%%%% Compute the y component of HEAT FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute using a different finite difference (choice will depend on which
% step in MacCormacks method we are in!)
% variables needed: ktherm, T
% case variable - needed so that can change which finite difference to use
% case_var should be a STRING!

% to compute along boundaries, use forward difference for left and lower,
% and backward difference for right and upper!!

% case 1: 'Predict_F' (predictor step computation)
% case 2: 'Correct_F' (corrector step computation)

function Qy = heatflux_y(delx,dely,ktherm,T,case_var)
Nx = size(T,2); Ny = size(T,1);
dTdy = zeros(size(T));

switch case_var
    %%%% Backward in y
    case 'Predict_F'
        %%%% INTERIOR AND BOTTOM BOUNDARY
        for i = 2:Ny
           
            dTdy(i,:) = (T(i,:)-T(i-1,:))/(dely);
            
        end         
        %%%% TOP BOUNDARY
        dTdy(1,:) = (T(2,:)-T(1,:))/(dely);
        
        Qy = -ktherm.*dTdy;
        
        
    %%%% Forward in y    
    case 'Correct_F'
        %%%% INTERIOR AND TOP BOUNDARY
        for i = 1:Ny-1
           
            dTdy(i,:) = (T(i+1,:)-T(i,:))/(dely);
            
        end         
        %%%% BOTTOM BOUNDARY
        dTdy(end,:) = (T(end,:)-T(end-1,:))/(dely);
        
        Qy = -ktherm.*dTdy;
        

%%%%%%%% Error case
    otherwise
        error('No finite difference scheme specified / incorrect "case_var"')
end
end
