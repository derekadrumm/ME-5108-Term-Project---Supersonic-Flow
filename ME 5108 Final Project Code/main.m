%%%% Main file for the final project %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project is to numerically solve the 2D Navier-Stokes equations for a
% supersonic flow over a flat plate using MacCormack Method
% EXTENDED DOMAIN!!
%
% We will solve for 9 unkowns:
% U,V, sqrt|U^2 + V^2|: velocity fields, velocity norm
% rho : density field
% P : pressure field
% T : temperature field
% Einternal : internal energy
% mu : viscosity
% ktherm : thermal conductivity
%
% **NOTE** Will need to add the 'function' folder to the file path!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%% Choose whether to extend domain past the plate
domain_case = 'Limited';
% domain_case = 'Extended';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch domain_case
case 'Extended'
% space grid
Nx = 141; % number of x points
Ny = 141; % number of y points

parameters; % problem parameters
% parameter vector which is used in many functions
params = [cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl];

left = 0.2*10^-5; right = left+Lplate; % ends of plate
up = 2; % extend domain in y
Ldomain = Lplate + 2*left;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otherwise
% space grid
Nx = 71; % number of x points
Ny = 71; % number of y points

parameters; % problem parameters
% parameter vector which is used in many functions
params = [cv, cp, R, mu0, ufar, Presfar, Tfar,Twall, prandtl];

left = 0; right = Lplate;
Ldomain = Lplate;
up = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delx = Ldomain/(Nx-1);
dely = (up*5*deltaBL)/(Ny-1); % y length needs to be 5*deltaBL to capture shock

x = 0:delx:delx*(Nx-1); y = 0:dely:dely*(Ny-1);
[X,Y] = meshgrid(x,y);

% Plate locations in the domain (used in boundary conditions)
global left_idx right_idx
left_idx = round(left/delx + 1);
right_idx = round(right/delx + 1);


initial_conditions; % intial U,V,P,T values, as well as rho, mu, ktherm


% Residual storage
Resid1 = [];
Resid2 = Resid1;
Resid3 = Resid1;
Resid4 = Resid1;

% method parameters
CFL = 0.6; % CFL number for time step computation
eps = 10E-12; % stop condition for the residual
iteration_stop = 8000; % stop condition for total iterations
residual = 1; iteration = 0; % initialized residual and iteration

%%%% MacCormack Method
tic
while residual > eps && iteration < iteration_stop
   
    iteration = iteration + 1;
    
    compute_timestep; % compute 'delt'
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Conservative Variable Computation (predictor step)
    [U1,U2,U3,U4] = conservative_U(rho,U,V,T,params);
    
    % When computing E and F, we need to use finite difference opposite to
    % what will be used in the predictor step computation!!
    
    % computations for E (backward x, central y)
    Qx = heatflux_x(delx,dely,ktherm,T,'Predict_E');
    
    Tauxx = stress_xx(delx,dely,U,V,mu,'Predict_E');
    Tauxy = stress_xy(delx,dely,U,V,mu,'Predict_E');
    
    [E1,E2,E3,E4] = conservative_E(rho,U,V,P,T,Tauxx,Tauxy,Qx,params);
    
    % computations for F (central x, backward y)
    Qy = heatflux_y(delx,dely,ktherm,T,'Predict_F');
    
    Tauxy = stress_xy(delx,dely,U,V,mu,'Predict_F');
    Tauyy = stress_yy(delx,dely,U,V,mu,'Predict_F');    
    
    [F1,F2,F3,F4] = conservative_F(rho,U,V,P,T,Tauxy,Tauyy,Qy,params);
    
    
    %%%% Predictor Step Computation (forward in space)
    U1pred = predictor_computation(delx,dely,delt,U1,E1,F1);
    U2pred = predictor_computation(delx,dely,delt,U2,E2,F2);
    U3pred = predictor_computation(delx,dely,delt,U3,E3,F3);
    U4pred = predictor_computation(delx,dely,delt,U4,E4,F4);
    
    
    %%%% Primitive Computation from Predictor
    [U,V,P,T] = primitive_UVPT(U1pred,U2pred,U3pred,U4pred,params);
    
    % Interpolated Boundary conditions
    [U,V,P,T] = boundary_conditions(U,V,P,T,params,domain_case);
    
    [mu,ktherm,rho] = prim_visc_therm_rho(P,T,params);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Conservative variable computations (Corrector step)
    % computations for E (forward x, central y)
    Qx = heatflux_x(delx,dely,ktherm,T,'Correct_E');
    
    Tauxx = stress_xx(delx,dely,U,V,mu,'Correct_E');
    Tauxy = stress_xy(delx,dely,U,V,mu,'Correct_E');
    
    [E1,E2,E3,E4] = conservative_E(rho,U,V,P,T,Tauxx,Tauxy,Qx,params);
    
    % computations for F (central x, forward y)
    Qy = heatflux_y(delx,dely,ktherm,T,'Correct_F');
    
    Tauxy = stress_xy(delx,dely,U,V,mu,'Correct_F');
    Tauyy = stress_yy(delx,dely,U,V,mu,'Correct_F');    
    
    [F1,F2,F3,F4] = conservative_F(rho,U,V,P,T,Tauxy,Tauyy,Qy,params);
    
    
    %%%% Corrector step computation (backward in space)
    U1cor = corrector_computation(delx,dely,delt,U1,U1pred,E1,F1);
    U2cor = corrector_computation(delx,dely,delt,U2,U2pred,E2,F2);
    U3cor = corrector_computation(delx,dely,delt,U3,U3pred,E3,F3);
    U4cor = corrector_computation(delx,dely,delt,U4,U4pred,E4,F4);
    
    
    %%%% Primitive Variable Solution
    [U,V,P,T] = primitive_UVPT(U1cor,U2cor,U3cor,U4cor,params);

    % Interpolated Boundary conditions
    [U,V,P,T] = boundary_conditions(U,V,P,T,params,domain_case);
    
    [mu,ktherm,rho] = prim_visc_therm_rho(P,T,params);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Residual computations
    Resid1 = [Resid1 norm(U1-U1cor,'fro')];
    Resid2 = [Resid2 norm(U2-U2cor,'fro')];
    Resid3 = [Resid3 norm(U3-U3cor,'fro')];
    Resid4 = [Resid4 norm(U4-U4cor,'fro')];
    
    residual = max([Resid1(end) Resid2(end) Resid3(end) Resid4(end)]);

    [Rex,Rey] = cell_reynolds(delx,dely,U,V,rho,mu);

end
comp_time = toc;

if iteration >= iteration_stop
    disp('Kill time reached')
elseif residual <= eps
    disp('residual epsilon reached')
end
disp(['MacCormack total computation time: ',num2str(comp_time)])


%%%% Plots and Error Analysis
yidx = length(y(y<5*deltaBL)); % index used for determining plot size

% Temperature Contour
figure(1)
hlabs = [300 340 420];
[C,h] = contourf(X(1:yidx,:),Y(1:yidx,:),T(1:yidx,:));
% [C,h] = contourf(X,Y,T);
clabel(C,h,hlabs,'Color','k')
colorbar
title('Temperature Field')

% Pressure Contour
figure(2)
contourf(X(1:yidx,:),Y(1:yidx,:),P(1:yidx,:));
colorbar
title('Pressure Field')

% Density Contour
figure(3)
contourf(X(2:yidx,:),Y(2:yidx,:),rho(2:yidx,:)); % for shock visualization
% contourf(X,Y,rho);
colorbar
title('Density Field')


% Residual semilogy
figure(4)
semilogy(1:iteration,Resid1,1:iteration,Resid2,1:iteration,Resid3,1:iteration,Resid4)
legend('U1','U2','U3','U4')
title('Semilogy Plot of Residuals of Conservative Variables')


% Shear Stress and Heat Transfer along plate
shear_stress_plate;
heat_transfer_plate;

drag_plate = trapz(tau_shear(left_idx:right_idx))
heatrate_plate = trapz(heat_plate(left_idx:right_idx))

figure(5)
plot(x,tau_shear)
xlabel('Length of Plate')
title('Shear Stress')

figure(6)
plot(x,heat_plate)
xlabel('Length of Plate')
title('Heat Transfer')


% Mach Number Field
% mach number m given by velocity/sqrt(gamma*R*T)
mach_field = sqrt(U.^2 + V.^2)./sqrt(gamma*R.*T);
mach_field = round(mach_field,5);

figure(7)
contourf(X(1:yidx,:),Y(1:yidx,:),mach_field(1:yidx,:))
colorbar
title('Mach Number Field','Rounded to nearest 10^{-4}')


