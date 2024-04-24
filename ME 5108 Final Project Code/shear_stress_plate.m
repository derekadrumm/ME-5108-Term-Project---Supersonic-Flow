%%%% SHEAR STRESS ALONG PLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% goal to analyze shear stress along the length of the plate
% shear stress along plate given by:
% 
% tau_shear = mu*(du/dy) evaluated along y=0
%
% Using a second order accurate forward difference:
%
% df/dx = (1/(2*delx))*(-3f_i + 4f_i+1 - f_i+2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% du/dy at y=0
dudy0 = (1/2*dely)*(-3*U(1,:) + 4*U(2,:) - U(3,:));

tau_shear = mu(1,:).*dudy0;
