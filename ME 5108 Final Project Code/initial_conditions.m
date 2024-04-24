%%%% Initial and Boundary Condition File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: some of these lines are unnecessary, but they are there to
% illistrate the initial/boundary conditions explicitly

ONE = ones(size(X)); % matrix of ones for initialization

T = ONE*Tfar;
P = ONE*Presfar;
U = ONE*ufar;
V = ONE*vfar;

[U,V,P,T] = boundary_conditions(U,V,P,T,params,domain_case);


% Initial Local Viscocity ()
mu = (mu0*((T/Tfar).^(3/2))).*((Tfar+110)./(T+110)); % local viscocity (need to compute this after each iteration!)

% Initial Local Density ()
rho = P./(R*T); % local density (need to compute this after each iteration!)

% Initial internal energy ()
Einternal = cv*T; % local internal energy

% Thermal Conductivity ()
ktherm = (mu*cp)/prandtl;
