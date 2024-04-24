%%%% COMPUTE NEW BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of new boundary conditions, which change for the outflow
% and plate surface at each time step!
% different cases for each boundary:
% 
% case_var: denotes if in limited or extended domain

function [U,V,P,T] = boundary_conditions(U,V,P,T,params,case_var)
global left_idx right_idx
paramcell = num2cell(params);
[cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl] = paramcell{:};

switch case_var
%%%% IF IN THE LIMITED DOMAIN    
case 'Limited'
% Case 1: at left corner of plate (i=1,j=1)
T(1,1) = Tfar;
P(1,1) = Presfar;
U(1,1) = 0;
V(1,1) = 0;

% Case 2: inflow and upper boundary (i=2:end,j=1) / (i=end,j=:)
T(2:end,1) = Tfar;      T(end,:) = Tfar;
P(2:end,1) = Presfar;   P(end,:) = Presfar;
U(2:end,1) = ufar;      U(end,:) = ufar;
V(2:end,1) = 0;         V(end,:) = 0;

% Case 3: plate surface (i=1,j=2:end)
T(1,2:end) = Twall; % temperature of plate
U(1,2:end) = 0; % no slip
V(1,2:end) = 0; % no slip
P(1,2:end) = 2*P(2,2:end) - P(3,2:end); % extrapolated from above plate

% Case 4: outflow, not surface or upper boundary (i=end,j=2:end-1)
T(2:end-1,end) = 2*T(2:end-1,end-1) - T(2:end-1,end-2); % extrapolated from interior
P(2:end-1,end) = 2*P(2:end-1,end-1) - P(2:end-1,end-2); % "
U(2:end-1,end) = 2*U(2:end-1,end-1) - U(2:end-1,end-2); % "
V(2:end-1,end) = 2*V(2:end-1,end-1) - V(2:end-1,end-2); % "

case 'Extended'
%%%% IF IN THE EXTENDED DOMAIN
% in this case, we use a symmetry boundary condition along y=0, where we 
% assume the gradients at both boundaries are 0. We compute the gradients 
% using a second order finite difference

% Case 1: at left corner of plate (i=1,j=1)
T(1,left_idx) = Tfar;
P(1,left_idx) = Presfar;
U(1,left_idx) = 0;
V(1,left_idx) = 0;

% Case 2: inflow and upper boundary (i=2:end,j=1) / (i=end,j=:)
T(2:end,1) = Tfar;      T(end,:) = Tfar;
P(2:end,1) = Presfar;   P(end,:) = Presfar;
U(2:end,1) = ufar;      U(end,:) = ufar;
V(2:end,1) = 0;         V(end,:) = 0;

% Case 3: plate surface (i=1,j=2:end)
T(1,left_idx:right_idx) = Twall; % temperature of plate
U(1,left_idx:right_idx) = 0; % no slip
V(1,left_idx:right_idx) = 0; % no slip
P(1,left_idx:right_idx) = 2*P(2,left_idx:right_idx) - P(3,left_idx:right_idx); % extrapolated from above plate

% Case 4: outflow (i=end,j=2:end-1)
T(2:end-1,end) = 2*T(2:end-1,end-1) - T(2:end-1,end-2); % extrapolated from interior
P(2:end-1,end) = 2*P(2:end-1,end-1) - P(2:end-1,end-2); % "
U(2:end-1,end) = 2*U(2:end-1,end-1) - U(2:end-1,end-2); % "
V(2:end-1,end) = 2*V(2:end-1,end-1) - V(2:end-1,end-2); % "

% Case 5: symmetry condition at y=0 (not plate surface)
% d/dy = 0 along y=0, so using 3 point finite difference
% ==> f_i = (4f_i+1 - f_i+2)/3

T(1,[1:left_idx-1,right_idx+1:end]) = (1/3)*(4*T(2,[1:left_idx-1,right_idx+1:end]) - T(3,[1:left_idx-1,right_idx+1:end]));
P(1,[1:left_idx-1,right_idx+1:end]) = (1/3)*(4*P(2,[1:left_idx-1,right_idx+1:end]) - P(3,[1:left_idx-1,right_idx+1:end]));
U(1,[1:left_idx-1,right_idx+1:end]) = (1/3)*(4*U(2,[1:left_idx-1,right_idx+1:end]) - U(3,[1:left_idx-1,right_idx+1:end]));
V(1,[1:left_idx-1,right_idx+1:end]) = 0; % velocity normal to bottom boundary is 0!


%%%%%%%% Error case
otherwise
    error('No finite difference scheme specified / incorrect "case_var"')
end






end
