%%%% Corrector step computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the INTERIOR of corrected UIcor
% given variables: UI,UIpred,EIpred,FIpred,delx,dely,delt
% partial derivatives in space computed using BACKWARD DIFFERENCE

function UIcor = corrector_computation(delx,dely,delt,UI,UIpred,EIpred,FIpred)
Nx = size(UI,2); Ny = size(UI,1);
UIcor = UIpred;

%%%% INTERIOR
for i = 2:Ny-1
    for j = 2:Nx-1
                
        UIcor(i,j) = (1/2)*(UI(i,j) + UIpred(i,j) - ...
            (delt/delx)*(EIpred(i,j)-EIpred(i,j-1)) - ...
            (delt/dely)*(FIpred(i,j)-FIpred(i-1,j)));
                
    end
end
end