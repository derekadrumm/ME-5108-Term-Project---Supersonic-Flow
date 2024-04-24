%%%% Parameter definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Fixed Problem Parameters (shouldnt change these)
% parameter = value; % Parameter name/definition (units)
Lplate = 0.00001; % length of the plate (m)
prandtl = 0.71; % Prandtl number for calorically perfect air ()
mu0 = 1.7894E-5; % reference viscocity at sea level (kg/(m s))
Presfar = 101325; % free stream pressure at sea level ()
Tfar = 288.15; % reference temperature at sea level (K)
R = 287; % perfect gas constant (J/(kg K))
gamma = 1.4; % ()
Reynolds = 1000; % Reynolds number (dimensionless)
Courant = 0.6; % Courant number should be 0.5 <= C <= 0.8 (dimensionless)


%%%% Variable Problem Parameters
Twall = 288.15; % temperature of plate
Machfar = 10; % desired mach number



%%%% Parameter Calculations (shouldnt change these)
cv = R/(gamma-1); % specific volume of air at sea level ()
cp = gamma*cv; % specific heat of air at sea level ()

soundspd = sqrt(gamma*R*Tfar); % speed of sound (m/s)

ufar = Machfar*soundspd; % free stream velocity (m/s)
vfar = 0; % farfield y velocity (0 by assumption) (m/s)
rhofar = Presfar/(R*Tfar); % density of air at sea level (kg/m^3)

Reynolds = (rhofar*ufar*Lplate)/mu0;
deltaBL = 5*Lplate/sqrt(Reynolds); % boundary layer thickness

