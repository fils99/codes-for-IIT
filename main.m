% se calcoli errore, commenta le prime 4 righe

% close all
% clear all
% clc
 
   area_max = 0.002; % attiva questo quando calcoli errore temporale 
% istanti_valori = 40; % attiva questo quando calcoli errore spaziale 

global epsilon geom  u_d  b_tot beta gamma grad_phi U_gdl u 
global phi_1 phi_2 phi_3 A_D_tot A_tot_lumping Pe_h tau_h R_tot_lumping
global istanti T_fin t delta_t Ndof M M_D term_noto b_N


istanti = istanti_valori;
T_fin = 0.5;
t = linspace(0, T_fin, istanti);
delta_t = t(2) - t(1);

Sample_Square_Dirichlet 

% assicurarsi che questi coefficienti rispettino la coercività della forma
% bilineare
epsilon = @ (x,y) 10^(0); 
beta = @ (x,y) [10*x;10*y] ;
gamma = @ (x,y) 20*x + y;

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
u = zeros (Np,istanti);


U_gdl = zeros(Ndof,istanti);

for j = 1: Np
    jj = geom.pivot.pivot(j);
    if jj > 0
        U_gdl(jj,1) = cond_iniz(geom.elements.coordinates(j,1), geom.elements.coordinates(j,2));
    end
end


% U_gdl = A_tot_lamping \ (b_tot - A_D_tot * u_d');

for n = 1:(istanti-1)
    U_gdl(:,n+1) = (M + delta_t/2 * A_tot_lumping ) \ (M*U_gdl(:,n) - M_D * (u_d(:,n+1) - u_d(:,n) )...
                 - delta_t/2 * A_tot_lumping * U_gdl(:,n) - delta_t/2 * A_D_tot * (u_d(:,n+1) + u_d(:,n))...
                 + delta_t /2 * (b_tot(:,n) + b_tot(:,n+1) )); 
end

for j = 1 : Np
    jj = geom.pivot.pivot(j);
    if jj > 0
        u (j,:) = U_gdl(jj,:);
    elseif jj<0
        u(j,:) = u_d(-jj,:);
    end
end

u_ex = zeros(Np,istanti);
for n = 1:istanti
u_ex(:,n) = sol_esatta_bis (geom.elements.coordinates(:,1), geom.elements.coordinates(:,2),t(n)); % sol esatta o sol esatta bis
end


if r < lung
    close(figure(1))
end

if r == lung

    % plot soluzione approssimata

    for n=1:istanti
    figure(2)
    plotta_soluz(geom.elements.triangles, geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u(:,n) )
    hold on
    title ('soluzione approssimata istante ',n)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off
    pause(0.01)
    end

    % plot soluzione esatta

    for n=1:istanti
    figure(3)
    plotta_soluz(geom.elements.triangles, geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u_ex(:,n) )
    hold on
    title ('soluzione esatta istante ',n)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off
    pause(0.01)
    end

end



%  L_Inf =  max ( abs ( u - u_ex ) );
%  
% 
% max(L_Inf)
% min(L_Inf)


% figure(2)
%     plotta_soluz(geom.elements.triangles, geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u(:,n) )
%     hold on
%     title ('soluzione approssimata istante ',n)
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     hold off
%     pause(0.01)
%     
% 
%     % plot soluzione esatta
% 
%     
%     figure(3)
%     plotta_soluz(geom.elements.triangles, geom.elements.coordinates(:,1), geom.elements.coordinates(:,2), u_ex(:,n) )
%     hold on
%     title ('soluzione esatta istante ',n)
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     hold off


