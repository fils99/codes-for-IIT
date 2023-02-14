% se calcoli errore, commenta le prime 4 righe

% close all
% clear all
% clc
% 
% area_max = 0.02;

global epsilon geom  u_d  b_tot beta gamma grad_phi U_gdl u phi_1 phi_2 phi_3 A_D_tot A_tot_lumping Pe_h tau_h R_tot_lumping


Sample_Square_Dirichlet 

% assicurarsi che questi coefficienti rispettino la coercività della forma
% bilineare
epsilon = @ (x,y) 10^(0); 
beta = @ (x,y) [0;0] ;
gamma = @ (x,y) 10000;

% funzioni test e loro gradienti
phi_1 = @(x,y) x;
phi_2 = @(x,y) y; 
phi_3 = @(x,y) (1-x-y);

grad_phi_1 = [1;0];
grad_phi_2 = [0;1];
grad_phi_3 = [-1;-1];

grad_phi = [grad_phi_1, grad_phi_2, grad_phi_3];

assembla_matrice ()




Np = geom.nelements.nVertexes;
u = zeros (Np,1);
%x = pcg (A, b - A_D * u_d');
 U_gdl = A_tot_lumping \ (b_tot - A_D_tot * u_d');

for j = 1 : Np
    jj = geom.pivot.pivot(j);
    if jj > 0
        u (j) = U_gdl(jj);
    elseif jj<0
        u(j) = u_d(-jj);
    end
end

u_ex = sol_esatta_bis (geom.elements.coordinates(:,1), geom.elements.coordinates(:,2)); % sol_esatta o sol_esatta_bis


if r < lung
    close(figure(1))
end

if r == lung

    % plot soluzione approssimata
    figure(2)
    
    plotta_soluz(geom.elements.triangles, geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u )
    hold on
    title ('soluzione approssimata')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    hold off

    % plot soluzione esatta
   
     figure(3)
    plotta_soluz(geom.elements.triangles, geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u_ex )
    hold on
    title ('soluzione esatta')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    hold off

end



 L_Inf =  max ( abs ( u - u_ex ) );





