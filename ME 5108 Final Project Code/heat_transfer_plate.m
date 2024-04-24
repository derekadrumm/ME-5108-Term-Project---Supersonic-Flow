%%%% HEAT TRANSFER ALONG PLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% goal to analyze heat transfer along the length of the plate
% heat transfer along plate given by:
% 
% heat_plate = ktherm*(dT/dy) evaluated along y=0
%
% Using a second order accurate forward difference:
%
% df/dx = (1/(2*delx))*(-3f_i + 4f_i+1 - f_i+2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dT/dy at y=0
dTdy0 = (1/2*dely)*(-3*T(1,:) + 4*T(2,:) - T(3,:));

heat_plate = ktherm(1,:).*dTdy0;
