%%%% Compute the x component of HEAT FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute using a different finite difference (choice will depend on which
% step in MacCormacks method we are in!)
% variables needed: ktherm, T
% case variable - needed so that can change which finite difference to use
% case_var should be a STRING!

% to compute along boundaries, use forward difference for left and lower,
% and backward difference for right and upper!!

% case 1: 'Predict_E' (predictor step computation)
% case 2: 'Correct_E' (corrector step computation)

function Qx = heatflux_x(delx,dely,ktherm,T,case_var)
Nx = size(T,2); Ny = size(T,1);
dTdx = zeros(size(T));

switch case_var
    %%%% Backward in x
    case 'Predict_E'
        %%%% INTERIOR AND RIGHT BOUNDARY
        for j = 2:Nx
           
            dTdx(:,j) = (T(:,j)-T(:,j-1))/(delx);
            
        end         
        %%%% LEFT BOUNDARY
        dTdx(:,1) = (T(:,2)-T(:,1))/(delx);
        
        Qx = -ktherm.*dTdx;
        
        
    %%%% Forward in x     
    case 'Correct_E'
        %%%% INTERIOR AND LEFT BOUNDARY
        for j = 1:Nx-1
           
            dTdx(:,j) = (T(:,j+1)-T(:,j))/(delx);
            
        end         
        %%%% RIGHT BOUNDARY
        dTdx(:,end) = (T(:,end)-T(:,end-1))/(delx);
        
        Qx = -ktherm.*dTdx;
        

%%%%%%%% Error case
    otherwise
        error('No finite difference scheme specified / incorrect "case_var"')
end
end
