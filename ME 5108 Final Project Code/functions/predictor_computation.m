%%%% Predictor step computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the INTERIOR of predicted UIpred
% given variables: UI,EI,FI,delx,dely,delt
% spacial derivatives computed using FORWARD DIFFERENCE

function UIpred = predictor_computation(delx,dely,delt,UI,EI,FI)
Nx = size(UI,2); Ny = size(UI,1);
% UIpred = zeros(size(UI));
UIpred = UI;

%%%% INTERIOR
for i = 2:Ny-1
    for j = 2:Nx-1
                
        UIpred(i,j) = UI(i,j) - (delt/delx)*(EI(i,j+1)-EI(i,j)) - (delt/dely)*(FI(i+1,j)-FI(i,j));
                
    end
end
end
