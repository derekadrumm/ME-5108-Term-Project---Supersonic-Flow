%%%% Compute the time step size 'delt' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltnew = 1;
nuijnew = 0;
% compute nuij
for i = 2:Ny-1
    for j = 2:Nx-1
        
        nuij = (4/3)*(mu(i,j)^2)*gamma;
        nuij = nuij/(prandtl*rho(i,j));
        
        nuijnew = max(nuij,nuijnew);
        
    end
end
nuij = nuijnew;

for i = 2:Nx-1
    for j = 2:Ny-1
        
        deltold = deltnew;
        a = sqrt(gamma*R*T(i,j)); % local sound speed
        
        deltnew = abs(U(i,j))/delx + abs(V(i,j))/dely + a*sqrt(1/(delx^2) + ...
            1/(dely^2)) + 2*nuij*(1/(delx^2) + 1/(dely^2));
        deltnew = CFL/deltnew;
        
        deltnew = min(deltnew,deltold);
        
    end
end

delt = deltnew;